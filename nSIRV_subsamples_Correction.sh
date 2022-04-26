#!/bin/bash

#RUN conda activate isoncorrect

#Author: Alexander Petri
# Modifications by Kristoffer Sahlin

# This shell program is used to run longterm tests on the algorithm. 
# For each amount of isoforms, this script performs 5 individual runs of 
# generating test data and running the algorithm on it. To find the remaining 
# bugs the script also outputs the test data files, which the algorithm did 
# not work on correctly

# RUN script as: ./generateTestResults.sh /path/to/input/reference.fa /path/to/output/folder

if [ $# -lt 1 ]; then
    # TODO: print usage
    echo " </path/to/input/reference.fa> <output_root>"
    exit 1
fi

input_ref=$1
filedirectory=$2
echo $filedirectory
mkdir -p $filedirectory
mkdir -p $filedirectory/errors/
mkdir -p $filedirectory/reads/
mkdir -p $filedirectory/isonform/
outputfile=$filedirectory/resultserror1.tsv
#if results.tsv already exists 
if [ -s $outputfile ]
then 
#delete all data from the file 
   > $outputfile
else
#add a file results.tsv to write into
   touch $outputfile
fi
#counter used to name error files correctly
errorcounter=0
#write nice header into the output file
#echo -e "Number of Isoforms \t Found isoforms by IsONform \n">>results.tsv
#ls
#define the file we want to use as indicator for our algos performance
file=$filedirectory/isonform/mapping.txt
echo -e "# reads \t #resulting isoforms  \t" >> $outputfile
#iterate over different numbers of isoforms
for ((i=1; i<=10; i++))
do
	#we want to have some double reads 
	n_reads=$(($i*100))
	#echo "Generating $i TestIsoforms" >>results.tsv
	#for each amount of isoforms we would like to run 10 tests to make sure our algo works stable
	for((j=1;j<=5;j++))
	do
		echo $i_$j
		outputs=0
		#python generateTestCases.py --ref $input_ref --sim_genome_len 1344 --nr_reads 20 --outfolder testout --coords 50 100 150 200 250 300 350 400 450 500 --probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms 2 --e True
		#run the test case generation script with the parameters needed
		
		number="${i}_${j}"
		############ COMMENT THE FOLLOWING TWO LINES FOR BUGFIXING ON IDENTICAL READ FILES ############
		###############################################################################################
		#python generateTestCases.py --ref $input_ref --sim_genome_len 1344 --nr_reads $n_reads --outfolder $filedirectory/reads --coords 50 100 150 200 250 300 350 400 450 500 --probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms $i --e True
		seqtk sample -s $j ERR3588903_1.fastq $n_reads >$filedirectory/reads/reads.fq
		mv $filedirectory/reads/reads.fq $filedirectory/reads_$number.fq
		###############################################################################################
		###############################################################################################

		isONcorrect --fastq $filedirectory/reads_$number.fq --outfolder $filedirectory/corrections
		mv $filedirectory/corrections/corrected_reads.fastq $filedirectory/reads_corr_$number.fq
		# cp $filedirectory/reads/reads.fq $filedirectory
		#we want to figure out how many reads were actually generated
		#read_amount=$(< $filedirectory/reads_$number.fq wc -l)
		#As fastq entries have 4 lines divide by 4 to get the actual number of reads
		#var=4
		#true_read_amount=$((read_amount / var))
		#run IsONform
		#if e=True
		
		python -m pyinstrument main.py --fastq $filedirectory/reads_corr_$number.fq --k 9 --w 20 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 5 --outfolder $filedirectory/isonform/
		#if e=False
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --delta_len 3 --outfolder out
		had_issue=$?
		zero=0
		echo $had_issue
		#python main.py --fastq $filedirectory/reads/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout 
		#we count the lines in mapping.txt, which is an output of IsONform. As we have a fasta file we have to divide this number by 2 to retreive the actual number of isoforms we generated
		res=$(< "$file" wc -l)
		result=$(($res/2))
		if [[ "$had_issue" != "$zero" ]]
		then
		errorcounter=$((errorcounter+1))
		result=-1
		cp $filedirectory/reads/reads.fq $filedirectory/errors/error_${errorcounter}_$number.fq
    fi
		#else
		#if [[ "$result" == "$true_read_amount" ]]
		#then
		#cp $filedirectory/reads/reads.fq $filedirectory/errors/strange_error_${errorcounter}_$number.fq

		#elif [[ "$result" -gt "$otherres" ]]
		#then
		#cp $filedirectory/reads/reads.fq $filedirectory/errors/other_error_${errorcounter}_$number.fq

		#elif [[ "$result" -lt "$i" ]]
		#then
		#cp $filedirectory/reads/reads.fq $filedirectory/errors/isoform_error_${errorcounter}_$number.fq
		#fi
		#fi
		#echo "$result"
		echo -e "$n_reads \t $result  \t" >> $outputfile


		#num2=2
		#otherres=$((result * num2))




		#if we have found an example for which IsONform does mess up, save it to debug later
		#if [[ "$i" != "$result" ]]
		#then
		#echo "Error found!"
		#cp testout/isoforms.fa error_${errorcounter}.fa
		#errorcounter=$((errorcounter+1))
		#echo "New file created" >>resultserror.tsv
		#fi
	done
	#echo -e"\n">>$outputfile
done


