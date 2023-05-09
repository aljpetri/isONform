#!/bin/bash
set -euo pipefail

# Download a small test dataset (D. melanogaster genome and some reads from the SRA)

mkdir -p tests/drosophila
cd tests/drosophila

if ! [[ -e ref.fastq ]]; then
  curl https://github.com/aljpetri/isONform_analysis/blob/main/test_data/SIRV_50_isos_cl0.fastq > ref.fastq.tmp
  mv ref.fastq.tmp ref.fastq
fi
#for r in 1 2; do
#    f=reads.${r}.fastq.gz
#    if ! [[ -e ${f} ]]; then
#      set +o pipefail
#      wget -nv -O - ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR605/006/SRR6055476/SRR6055476_${r}.fastq.gz | zcat | head -n 400000 | gzip > ${f}.tmp
#      set -o pipefail
#      mv ${f}.tmp ${f}
#    fi
#done