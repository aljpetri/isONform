#!/bin/bash
set -e
#if [ $# -ne 6 ]; then
#        echo "Usage: `basename $0`  <raw_reads.fq>  <outfolder>  <num_cores> <isONform_folder> <iso_abundance> <mode> "
#        exit 0
#fi
#the pipeline can be run in different modes:
# full: the full pipeline is run (pychopper, isONclust,isONcorrect,isONform)
# pacbio: for PacBio data runs isONclust and isONform
# analysis: special mode for the analysis pipelines: only isONclust,isONcorrect and isONform are run
# only_isonform: only isONform is run
raw_reads=$1
outfolder=$2
num_cores=$3
isONform_folder=$4
iso_abundance=$5
mode=$6
echo "Running: " `basename $0` $raw_reads $outfolder $num_cores $isONform_folder $iso_abundance $mode
isonform_folder=${isONform_folder::-1}
mkdir -p $outfolder
#conda init bash
#conda activate /proj/snic2022-6-31/nobackup/alexp/conda_envs/isonform
if [ $mode == "full" ]
then
echo
echo "Will run pychopper (cdna_classifier.py), isONclust, isONcorrect and isONform. Make sure you have these tools installed."
echo "For installation see: https://github.com/ksahlin/isONcorrect#installation and  https://github.com/aljpetri/isONform"
echo

echo
echo "Running pychopper"
echo

pychopper  $raw_reads $outfolder/full_length.fq -t $num_cores

echo
echo "Finished pychopper"
echo

fi

echo
echo "Running isONclust"
echo
if [ $mode != "only_isonform" ]
then
if [ $mode != "pacbio" ] && [ $mode != "analysis" ]
then
/usr/bin/time -v isONclust  --t $num_cores  --ont --fastq $outfolder/full_length.fq \
             --outfolder $outfolder/clustering
/usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $outfolder/full_length.fq --outfolder  $outfolder/clustering/fastq_files
elif [ $mode == "analysis" ]
then
/usr/bin/time -v  isONclust  --t $num_cores  --ont --fastq $raw_reads \
             --outfolder $outfolder/clustering --k 8 --w 9
/usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $raw_reads --outfolder  $outfolder/clustering/fastq_files
else
/usr/bin/time -v  isONclust  --t $num_cores  --isoseq  --fastq $raw_reads \
             --outfolder $outfolder/clustering
/usr/bin/time -v isONclust write_fastq --N $iso_abundance --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $raw_reads --outfolder  $outfolder/clustering/fastq_files
fi

echo
echo "Finished isONclust"
echo
#conda activate isON311

if [ $mode != "pacbio" ]
then
  echo
  echo "Running isONcorrect"
  echo

 /usr/bin/time -v python3.11 run_isoncorrect --t $num_cores  --fastq_folder $outfolder/clustering/fastq_files  --outfolder $outfolder/correction/

  echo
  echo "Finished isONcorrect"
  echo
fi
fi
echo
echo "Merging reads back to single file. Corrected reads per cluster are still stored in: " $outfolder/correction/
echo

echo
echo "Running isONform"
echo
#if [ $mode != "pacbio" ]
#then
#  python3.11 $isONform_folder/isONform_parallel.py --fastq_folder $outfolder/correction/ --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 15 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance --split_wrt_batches --merge_sub_isoforms_3  --merge_sub_isoforms_5 --delta_iso_len_3 30 --delta_iso_len_5 50 --slow
#else
/usr/bin/time -v  python3.11 $isONform_folder/isONform_parallel.py --t $num_cores --fastq_folder $outfolder/clustering/fastq_files --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 10 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance --split_wrt_batches --merge_sub_isoforms_3  --merge_sub_isoforms_5 --delta_iso_len_3 30 --delta_iso_len_5 50 --slow --clustered
#fi
echo
echo "Finished isONform"
echo

echo
echo "Finished with pipeline and wrote corrected reads to: " $outfolder
echo

