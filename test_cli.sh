#!/bin/bash
set -ex

# Run the executable with the demo BAM file
./target/release/pbr demo/test.bam

echo "CLI test passed successfully!"
