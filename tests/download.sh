#!/bin/bash
set -euo pipefail

# Download a small test dataset (D. melanogaster genome and some reads from the SRA)

mkdir -p tests/drosophila
cd tests/drosophila
mkdir -p 0
cd 0

#if ! [[ -e corrected_reads.fastq ]]; then
  curl -o corrected_reads.fastq https://github.com/aljpetri/isONform_analysis/blob/main/test_data/SIRV_50_isos_cl0.fastq
  #mv ref.fastq.tmp corrected_reads.fastq
  #echo "$PWD"
  echo "$(head -20 corrected_reads.fastq)"
  #echo cat my_file.txt

#fi
if ! [[ -e ref.fasta ]]; then
  curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz | gunzip > ref.fasta.tmp
  mv ref.fasta.tmp ref.fasta
  echo "$(head -20 ref.fasta)"
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