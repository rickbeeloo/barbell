#!/bin/bash
# Script to run flamegraph profiling on demuxer

# Default values
FASTQ="${1:-../tests/sample.fastq}"
KIT="${2:-SQK-RBK114-96}"
ITERATIONS="${3:-10}"

echo "Running flamegraph profiling..."
echo "FASTQ: $FASTQ"
echo "Kit: $KIT"
echo "Iterations: $ITERATIONS"
echo ""

# Check if cargo-flamegraph is installed
if ! command -v cargo-flamegraph &> /dev/null; then
    echo "cargo-flamegraph is not installed."
    echo "Install it with: cargo install flamegraph"
    exit 1
fi

# Run flamegraph
cargo flamegraph --bin profile_demux -- "$FASTQ" "$KIT" "$ITERATIONS"

echo ""
echo "Flamegraph generated: flamegraph.svg"
echo "Open it in your browser to view the profiling results."

