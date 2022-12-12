#!/bin/bash
set -e
if [ $# -ne 5 ]; then
        echo "Usage: `basename $0`  <raw_reads.fq>  <outfolder>  <num_cores> "
        exit 0
fi

raw_reads=$1
outfolder=$2
num_cores=$3
isONform_folder=$4
iso_abundance=$5
echo "Running: " `basename $0` $raw_reads $outfolder $num_cores $isONform_folder $iso_abundance
#isonform_folder=${isONform_folder::-1}
mkdir -p $outfolder
echo "ISONfolder "$isONform_folder
echo
echo "Will run pychopper (cdna_classifier.py), isONclust, isONcorrect and isONform. Make sure you have these tools installed."
echo "For installation see: https://github.com/ksahlin/isONcorrect#installation and  https://github.com/aljpetri/isONform"
echo

#echo
#echo "Running pychopper"
#echo

#pychopper  $raw_reads $outfolder/full_length.fq -t $num_cores

#echo
#echo "Finished pychopper"
#echo



echo
echo "Running isONclust"
echo

isONclust  --t $num_cores  --ont --fastq $raw_reads \
             --outfolder $outfolder/clustering  --k 8 --w 9 --min_shared 3 #old k=10 oldw=12
isONclust write_fastq --N 1 --clusters $outfolder/clustering/final_clusters.tsv \
                      --fastq $raw_reads --outfolder  $outfolder/clustering/fastq_files
echo
echo "Finished isONclust"
echo



echo
echo "Running isONcorrect"
echo
#run_isoncorrect --t $num_cores  --fastq_folder $raw_reads  --outfolder $outfolder/correction/
#the following line is without isONclust:
run_isoncorrect --t $num_cores  --fastq_folder $outfolder/clustering/fastq_files  --outfolder $outfolder/correction/

echo
echo "Finished isONcorrect"
echo


echo
echo "Running isONform"
echo
#python3.11 $isONform_folder/isONform_parallel.py --fastq_folder $outfolder/correction/ --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 5 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance --split_wrt_batches
#python3.11 $isONform_folder/main.py --fastq $outfolder/correction/*/*.fastq --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 5 --outfolder $outfolder/singleisoforms --iso_abundance $iso_abundance
python3.11 $isONform_folder/isONform_parallel.py --fastq_folder $outfolder/correction/ --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 5 --outfolder $outfolder/isoforms --iso_abundance $iso_abundance --split_wrt_batches --merge_sub_isoforms_3 --merge_sub_isoforms_5 --delta_iso_len_3 5 --delta_iso_len_5 5 --slow
echo
echo "Finished isONform"
echo
# OPTIONAL BELOW TO MERGE ALL CORRECTED READS INTO ONE FILE
#touch $outfolder/all_corrected_reads.fq
#OUTFILES=$outfolder"/correction/"*"/corrected_reads.fastq"
#for f in $OUTFILES
#do
#  echo $f
#  cat $f >> $outfolder/all_corrected_reads.fq
#done

echo
echo "Finished with pipeline and wrote corrected reads to: " $outfolder/isoforms
echo
