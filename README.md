# pbr

[![Rust](https://github.com/brentp/pbr/actions/workflows/rust.yml/badge.svg)](https://github.com/brentp/pbr/actions/workflows/rust.yml)
drunk on [perbase](https://github.com/sstadick/perbase) pileups and [lua](https://github.com/khvzak/mlua/) expressions.

## Usage

### Basic Usage

```bash
pbr <BAM_FILE> -f <FASTA_FILE> -t <THREADS> --mate-fix --max-depth <MAX_DEPTH> -k <FLANKING>
```

### Examples

Single Base Reference (`-k 0`)

```bash
pbr test/simple_test.bam -f test/simple_test.fa -t 4 --mate-fix --max-depth 500000 -k 0
```

Multi-Base Reference Window (`-k 1`)

```bash
pbr test/simple_test.bam -f test/simple_test.fa -t 4 --mate-fix --max-depth 500000 -k 1
```

Filter rows where the reference base is `'A'`:

```bash
pbr test/simple_test.bam -f test/simple_test.fa -t 4 --mate-fix --max-depth 500000 -k 0 --pile-expression "return pile.ref_base == 'A'"
```

Filter rows where the reference window is `'AAG'`:

```bash
pbr test/simple_test.bam -f test/simple_test.fa -t 4 --mate-fix --max-depth 500000 -k 1 --pile-expression "return pile.ref_base == 'AAG'"
```

Filter reads with a specific tag, e.g., `RG:Z:test`:

```bash
pbr test/simple_test.bam -f test/simple_test.fa -t 4 --mate-fix --max-depth 500000 -k 0 --pile-expression "return read:tag('RG') == 'test'"
```

## Lua Filtering

The `--expression` (`-e`) option allows you to filter reads using Lua expressions. The `read` object exposes the following fields:

- `qname`: Query name
- ....

The `--pile-expression` (`-p`) option allows you to filter rows using Lua expressions. The `pile` object exposes the following fields:

- `depth`: Read depth
- `a`, `c`, `g`, `t`, `n`: Counts of each base
- `ins`, `del`, `ref_skip`, `fail`: Counts of insertions, deletions, reference skips, and failed reads
- `pos`: 0-based position
- `ref_base`: Reference base or window (e.g., `'A'` for `-k 0`, `'AAG'` for `-k 1`)

## Output

The output is a tab-separated file with the following columns:

- `#chrom`: Chromosome name
- `pos0`: 0-based position
- `ref_base`: Reference base or window (e.g., `'A'` for `-k 0`, `'AAG'` for `-k 1`)
- `depth`: Read depth
- `a`, `c`, `g`, `t`, `n`: Counts of each base
