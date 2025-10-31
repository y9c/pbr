//! # ParGranges
//!
//! Iterates over chunked genomic regions in parallel.
use anyhow::{anyhow, Context, Result};
use bio::io::bed;
use crossbeam::channel::{bounded, Receiver};
use lazy_static::lazy_static;
use log::*;
use num_cpus;
use rayon::prelude::*;
use rust_htslib::{
    bam::{HeaderView, IndexedReader, Read},
    bcf::{Read as bcfRead, Reader},
};
use rust_lapper::{Interval, Lapper};
use serde::Serialize;
use std::{convert::TryInto, path::PathBuf, thread};

const BYTES_INA_GIGABYTE: usize = 1024 * 1024 * 1024;

/// A modifier to apply to the channel size formular that is (BYTES_INA_GIGABYTE * channel_size_modifier) * threads / size_of(R::P)
/// 0.15 roughly corresponds to 1_000_000 PileupPosition objects per thread with some wiggle room.
pub const CHANNEL_SIZE_MODIFIER: f64 = 0.15;

/// The ideal number of basepairs each worker will receive. Total bp in memory at one time = `threads` * `chunksize`
pub const CHUNKSIZE: u32 = 1_000_000;

lazy_static! {
    /// CHANNEL_SIZE_MODIFIER as a str
    pub static ref CHANNEL_SIZE_MODIFIER_STR: String = CHANNEL_SIZE_MODIFIER.to_string();

    /// CHUNKSIZE as a str
    pub static ref CHUNKSIZE_STR: String = CHUNKSIZE.to_string();
}

/// RegionProcessor defines the methods that must be implemented to process a region
pub trait RegionProcessor {
    /// A vector of P make up the output of [`process_region`] and
    /// are values associated with each position.
    ///
    /// [`process_region`]: #method.process_region
    type P: 'static + Send + Sync + Serialize;

    /// A function that takes the tid, start, and stop and returns something serializable.
    /// Note, a common use of this function will be a `fetch` -> `pileup`. The pileup must
    /// be bounds checked.
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P>;
}

/// ParGranges holds all the information and configuration needed to launch the
/// [`ParGranges::process`].
///
/// [`ParGranges::process`]: #method.process
#[derive(Debug)]
pub struct ParGranges<R: 'static + RegionProcessor + Send + Sync> {
    /// Path to an indexed BAM / CRAM file
    reads: PathBuf,
    /// Optional reference file for CRAM
    ref_fasta: Option<PathBuf>,
    /// Optional path to a BED file to restrict the regions iterated over
    regions_bed: Option<PathBuf>,
    /// Optional path to a BCF/VCF file to restrict the regions iterated over
    regions_bcf: Option<PathBuf>,
    /// If `regions_bed` and or `regions_bcf` is specified, and this is true, merge any overlapping regions to avoid duplicate output.
    merge_regions: bool,
    /// Number of threads this is allowed to use, uses all if None
    threads: usize,
    /// The ideal number of basepairs each worker will receive. Total bp in memory at one time = `threads` * `chunksize`
    chunksize: u32,
    /// A modifier to apply to the channel size formular that is (BYTES_INA_GIGABYTE * channel_size_modifier) * threads / size_of(R::P)
    channel_size_modifier: f64,
    /// The rayon threadpool to operate in
    pool: rayon::ThreadPool,
    /// The implementation of [RegionProcessor] that will be used to process regions
    processor: R,
}

impl<R: RegionProcessor + Send + Sync> ParGranges<R> {
    /// Create a ParIO object
    ///
    /// # Arguments
    ///
    /// * `reads`- path to an indexed BAM/CRAM
    /// * `ref_fasta`- path to an indexed reference file for CRAM
    /// * `regions_bed`- Optional BED file path restricting the regions to be examined
    /// * `regions_bcf`- Optional BCF/VCF file path restricting the regions to be examined
    /// * `merge_regions` - If `regions_bed` and or `regions_bcf` is specified, and this is true, merge any overlapping regions to avoid duplicate output.
    /// * `threads`- Optional threads to restrict the number of threads this process will use, defaults to all
    /// * `chunksize`- optional argument to change the default chunksize of 1_000_000. `chunksize` determines the number of bases
    ///                each worker will get to work on at one time.
    /// * `channel_size_modifier`- Optional argument to modify the default size ration of the channel that `R::P` is sent on.
    ///                 formula is: ((BYTES_INA_GIGABYTE * channel_size_modifier) * threads) / size_of(R::P)
    /// * `processor`- Something that implements [`RegionProcessor`](RegionProcessor)
    pub fn new(
        reads: PathBuf,
        ref_fasta: Option<PathBuf>,
        regions_bed: Option<PathBuf>,
        regions_bcf: Option<PathBuf>,
        merge_regions: bool,
        threads: Option<usize>,
        chunksize: Option<u32>,
        channel_size_modifier: Option<f64>,
        processor: R,
    ) -> Self {
        let threads = if let Some(threads) = threads {
            threads
        } else {
            num_cpus::get()
        };

        // Keep two around for main thread and thread running the pool
        let threads = std::cmp::max(threads.checked_sub(2).unwrap_or(0), 1);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap();

        info!("Using {} worker threads.", threads);
        Self {
            reads,
            ref_fasta,
            regions_bed,
            regions_bcf,
            merge_regions,
            threads,
            chunksize: chunksize.unwrap_or(CHUNKSIZE),
            channel_size_modifier: channel_size_modifier.unwrap_or(CHANNEL_SIZE_MODIFIER),
            pool,
            processor,
        }
    }

    /// Process each region.
    ///
    /// This method splits the sequences in the BAM/CRAM header into `chunksize` * `self.threads` regions (aka 'super chunks').
    /// It then queries that 'super chunk' against the intervals (either the BED file, or the whole genome broken up into `chunksize`
    /// regions). The results of that query are then processed by a pool of workers that apply `process_region` to reach interval to
    /// do perbase analysis on. The collected result for each region is then sent back over the returned `Receiver<R::P>` channel
    /// for the caller to use. The results will be returned in order according to the order of the intervals used to drive this method.
    ///
    /// While one 'super chunk' is being worked on by all workers, the last 'super chunks' results are being printed to either to
    /// a file or to STDOUT, in order.
    ///
    /// Note, a common use case of this will be to fetch a region and do a pileup. The bounds of bases being looked at should still be
    /// checked since a fetch will pull all reads that overlap the region in question.
    pub fn process(self) -> Result<Receiver<R::P>> {
        let channel_size: usize = ((BYTES_INA_GIGABYTE as f64 * self.channel_size_modifier).floor()
            as usize
            / std::mem::size_of::<R::P>())
            * self.threads;
        info!(
            "Creating channel of length {:?} (* 120 bytes to get mem)",
            channel_size
        );
        let (snd, rxv) = bounded(channel_size);
        thread::spawn(move || {
            self.pool.install(|| {
                info!("Reading from {:?}", self.reads);
                let mut reader = IndexedReader::from_path(&self.reads).expect("Indexed BAM/CRAM");
                // If passed add ref_fasta
                if let Some(ref_fasta) = &self.ref_fasta {
                    reader.set_reference(ref_fasta).expect("Set ref");
                }
                // Get a copy of the header
                let header = reader.header().to_owned();
                let header_names: Vec<String> = (0..header.target_count())
                    .map(|tid| std::str::from_utf8(header.tid2name(tid)).unwrap().to_string())
                    .collect();

                // Work out if we are restricted to a subset of sites
                let bed_intervals = if let Some(regions_bed) = &self.regions_bed {
                    Some(
                        Self::bed_to_intervals(&header, regions_bed, self.merge_regions)
                            .expect("Parsed BED to intervals"),
                    )
                } else {
                    None
                };
                let bcf_intervals = if let Some(regions_bcf) = &self.regions_bcf {
                    Some(
                        Self::bcf_to_intervals(&header, regions_bcf, self.merge_regions)
                            .expect("Parsed BCF/VCF to intervals"),
                    )
                } else {
                    None
                };
                let restricted_ivs = match (bed_intervals, bcf_intervals) {
                    (Some(bed_ivs), Some(bcf_ivs)) => {
                        Some(Self::merge_intervals(bed_ivs, bcf_ivs, self.merge_regions))
                    }
                    (Some(bed_ivs), None) => Some(bed_ivs),
                    (None, Some(bcf_ivs)) => Some(bcf_ivs),
                    (None, None) => None,
                };

                let intervals = if let Some(restricted) = restricted_ivs {
                    restricted
                } else {
                    Self::header_to_intervals(&header, self.chunksize)
                        .expect("Parsed BAM/CRAM header to intervals")
                };

                // Group small chromosomes together until reaching chunk size
                let mut current_batch = Vec::new();
                let mut current_batch_size = 0;

                for (tid, intervals) in intervals.into_iter().enumerate() {
                    let tid: u32 = tid as u32;
                    let tid_end: u32 = header.target_len(tid).unwrap().try_into().unwrap();
                    let chr_name = &header_names[tid as usize];
                    debug!("Checking chromosome: {} (TID: {})", chr_name, tid);
                    // Skip chromosomes with no reads
                    if !self.has_reads(&mut reader, tid) {
                        debug!("Skipping chromosome: {} (TID: {}) as it has no reads", chr_name, tid);
                        continue;
                    }

                    // If this is a large chromosome, process it in chunks
                    if tid_end > self.chunksize {
                        // Process any pending batch first
                        if !current_batch.is_empty() {
                            self.process_chromosome_batch(&mut reader, &header_names, &current_batch, &snd);
                            current_batch.clear();
                            current_batch_size = 0;
                        }

                        info!("Processing chromosome {}:0-{} (chunked)", tid, tid_end);
                        let mut result = vec![];
                        for chunk_start in (0..tid_end).step_by(self.chunksize as usize) {
                            let chunk_end = std::cmp::min(chunk_start as u32 + self.chunksize, tid_end);
                            let tid_name = &header_names[tid as usize];
                            trace!(
                                "Batch Processing {}:{}-{}",
                                tid_name,
                                chunk_start,
                                chunk_end
                            );

                            let (r, _) = rayon::join(
                                || {
                                    let ivs: Vec<Interval<u32, ()>> =
                                        Lapper::<u32, ()>::find(&intervals, chunk_start, chunk_end)
                                            .map(|iv| Interval {
                                                start: std::cmp::max(iv.start, chunk_start),
                                                stop: std::cmp::min(iv.stop, chunk_end),
                                                val: (),
                                            })
                                            .collect();
                                    ivs.into_par_iter()
                                        .flat_map(|iv| {
                                            trace!("Processing {}:{}-{}", tid_name, iv.start, iv.stop);
                                            self.processor.process_region(tid, iv.start, iv.stop)
                                        })
                                        .collect()
                                },
                                || {
                                    result.into_iter().for_each(|p| {
                                        snd.send(p).expect("Sent a serializable to writer")
                                    })
                                },
                            );
                            result = r;
                        }
                        result.into_iter().for_each(|p| snd.send(p).expect("Sent a serializable to writer"));
                    } else {
                        // This is a small chromosome, add it to the current batch
                        current_batch.push((tid, tid_end, intervals.clone()));
                        current_batch_size += tid_end;

                        // Process the batch if we've reached chunk size or finished all chromosomes
                        if current_batch_size >= self.chunksize || tid == intervals.len() as u32 - 1 {
                            self.process_chromosome_batch(&mut reader, &header_names, &current_batch, &snd);
                            current_batch.clear();
                            current_batch_size = 0;
                        }
                    }
                }

                // Process any remaining batch after going through all chromosomes
                if !current_batch.is_empty() {
                    self.process_chromosome_batch(&mut reader, &header_names, &current_batch, &snd);
                }
            });
        });
        Ok(rxv)
    }

    /// Process a batch of small chromosomes together
    fn process_chromosome_batch(
        &self,
        _reader: &mut IndexedReader,
        header_names: &[String],
        batch: &[(u32, u32, Lapper<u32, ()>)],
        snd: &crossbeam::channel::Sender<R::P>,
    ) {
        if batch.is_empty() {
            return;
        }

        let mut result = vec![];
        let (r, _) = rayon::join(
            || {
                batch.into_par_iter()
                    .flat_map(|(tid, tid_end, intervals)| {
                        let tid_name = &header_names[*tid as usize];
                        info!("Processing chromosome {}:0-{}", tid, tid_end);
                        let ivs: Vec<Interval<u32, ()>> = Lapper::<u32, ()>::find(intervals, 0, *tid_end)
                            .map(|iv| Interval {
                                start: iv.start,
                                stop: iv.stop,
                                val: (),
                            })
                            .collect();
                        ivs.into_par_iter()
                            .flat_map(|iv| {
                                trace!("Processing {}:{}-{}", tid_name, iv.start, iv.stop);
                                self.processor.process_region(*tid, iv.start, iv.stop)
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect()
            },
            || {
                result.into_iter().for_each(|p| {
                    snd.send(p).expect("Sent a serializable to writer")
                })
            },
        );
        result = r;
        result.into_iter().for_each(|p| snd.send(p).expect("Sent a serializable to writer"));
    }

    // Convert the header into intervals of equally sized chunks. The last interval may be short.
    fn header_to_intervals(header: &HeaderView, chunksize: u32) -> Result<Vec<Lapper<u32, ()>>> {
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for tid in 0..(header.target_count()) {
            let tid_len: u32 = header.target_len(tid).unwrap().try_into().unwrap();
            for start in (0..tid_len).step_by(chunksize as usize) {
                let stop = std::cmp::min(start as u32 + chunksize, tid_len);
                intervals[tid as usize].push(Interval {
                    start: start as u32,
                    stop: stop,
                    val: (),
                });
            }
        }
        Ok(intervals.into_iter().map(|ivs| Lapper::new(ivs)).collect())
    }

    /// Read a bed file into a vector of lappers with the index representing the TID
    /// if `merge' is true then any overlapping intervals in the sets will be merged.
    // TODO add a proper error message
    fn bed_to_intervals(
        header: &HeaderView,
        bed_file: &PathBuf,
        merge: bool,
    ) -> Result<Vec<Lapper<u32, ()>>> {
        let mut bed_reader = bed::Reader::from_file(bed_file)?;
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for (i, record) in bed_reader.records().enumerate() {
            let record = record?;
            let tid = header
                .tid(record.chrom().as_bytes())
                .expect("Chromosome not found in BAM/CRAM header");
            let start = record
                .start()
                .try_into()
                .with_context(|| format!("BED record {} is invalid: unable to parse start", i))?;
            let stop = record
                .end()
                .try_into()
                .with_context(|| format!("BED record {} is invalid: unable to parse stop", i))?;
            if stop < start {
                return Err(anyhow!("BED record {} is invalid: stop < start", i));
            }
            intervals[tid as usize].push(Interval {
                start,
                stop,
                val: (),
            });
        }

        Ok(intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect())
    }

    /// Read a BCF/VCF file into a vector of lappers with index representing the TID
    /// if `merge' is true then any overlapping intervals in the sets will be merged.
    fn bcf_to_intervals(
        header: &HeaderView,
        bcf_file: &PathBuf,
        merge: bool,
    ) -> Result<Vec<Lapper<u32, ()>>> {
        let mut bcf_reader = Reader::from_path(bcf_file).expect("Error opening BCF/VCF file.");
        let bcf_header_reader = Reader::from_path(bcf_file).expect("Error opening BCF/VCF file.");
        let bcf_header = bcf_header_reader.header();
        let mut intervals = vec![vec![]; header.target_count() as usize];
        // TODO: validate the headers against eachother
        for record in bcf_reader.records() {
            let record = record?;
            let record_rid = bcf_header.rid2name(record.rid().unwrap()).unwrap();
            let tid = header
                .tid(record_rid)
                .expect("Chromosome not found in BAM/CRAM header");
            let pos: u32 = record
                .pos()
                .try_into()
                .expect("Got a negative value for pos");
            intervals[tid as usize].push(Interval {
                start: pos,
                stop: pos + 1,
                val: (),
            });
        }

        Ok(intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect())
    }

    /// Merge two sets of restriction intervals together
    /// if `merge' is true then any overlapping intervals in the sets will be merged.
    fn merge_intervals(
        a_ivs: Vec<Lapper<u32, ()>>,
        b_ivs: Vec<Lapper<u32, ()>>,
        merge: bool,
    ) -> Vec<Lapper<u32, ()>> {
        let mut intervals = vec![vec![]; a_ivs.len()];
        for (i, (a_lapper, b_lapper)) in a_ivs.into_iter().zip(b_ivs.into_iter()).enumerate() {
            intervals[i] = a_lapper.into_iter().chain(b_lapper.into_iter()).collect();
        }
        intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect()
    }

    /// Check if a chromosome has any reads by fetching the entire chromosome
    fn has_reads(&self, reader: &mut IndexedReader, tid: u32) -> bool {
        reader.fetch(tid).expect("Failed to fetch chromosome");
        reader.pileup().next().is_some()
    }
}
