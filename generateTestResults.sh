#!/bin/bash
#Author: Alexander Petri
#This shell program is used to run longterm tests on the algorithm. For each amount of isoforms, this script performs 5 individual runs of generating test data and running the algorithm on it. To find the remaining bugs the script also outputs the test data files, which the algorithm did not work on correctly
mkdir -p filedirectory
outputfile=resultserror7.tsv
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
file=out/mapping.txt

#iterate over different numbers of isoforms
for ((i=3; i<=25; i++))
do
#we want to have some double reads 
n_reads=$(($i+5))
#echo "Generating $i TestIsoforms" >>results.tsv
#for each amount of isoforms we would like to run 5 tests to make sure our algo works stable
for((j=1;j<=5;j++))
do
echo $i_$j
outputs=0
#run the test case generation script with the parameters needed
python generateTestCases.py --ref /home/alexanderpetri/Desktop/RAWDATA_PhD1/Isoform_Test_data.fa --sim_genome_len 1344 --nr_reads $n_reads --outfolder testout --coords 50 100 150 200 250 300 350 400 450 500 --probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms $i --e True
cp ~/PHDProject1/testout/reads.fq filedirectory
#we want to figure out how many reads were actually generated
read_amount=$(< ~/PHDProject1/testout/reads.fq wc -l)
#As fastq entries have 4 lines divide by 4 to get the actual number of reads
var=4
true_read_amount=$((read_amount / var))
number="${i}_${j}"
echo $number
mv ~/PHDProject1/filedirectory/reads.fq ~/PHDProject1/filedirectory/reads_$number.fq
#run IsONform
python main.py --fastq ~/PHDProject1/testout/reads.fq --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout 
#python main.py --fastq ~/PHDProject1/testout/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout 
#we count the lines in mapping.txt, which is an output of IsONform. As we have a fasta file we have to divide this number by 2 to retreive the actual number of isoforms we generated
res=$(< "$file" wc -l)
result=$(($res/2))
#echo "$result"
echo -e "$i \t $result \t $true_read_amount \t" >> $outputfile
#if we have found an example for which IsONform does mess up, save it to debug later
#if [[ "$i" != "$result" ]]
#then
#echo "Error found!"
#cp testout/isoforms.fa error_${errorcounter}.fa
#errorcounter=$((errorcounter+1))
#echo "New file created" >>resultserror.tsv
#fi
done
echo -e"\n">>$outputfile
done

