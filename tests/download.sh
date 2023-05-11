#!/bin/bash
set -euo pipefail

# Download a small test dataset (D. melanogaster genome and some reads from the SRA)

mkdir -p tests/drosophila
cd tests/drosophila
mkdir -p 0
cd 0

if ! [[ -e corrected_reads.fastq ]]; then
  curl -o corrected_reads.fastq https://raw.githubusercontent.com/aljpetri/isONform_analysis/main/test_data/SIRV_50_isos_cl0.fastq
fi