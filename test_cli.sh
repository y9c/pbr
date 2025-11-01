#!/bin/bash
set -ex

# Run the executable with the demo BAM file
./target/release/pbr demo/test.bam -e "return read.mapq > 10"

echo "CLI test passed successfully!"
