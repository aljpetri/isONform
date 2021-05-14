#!/bin/bash
if [ -s results.tsv ]
then 
   > results.tsv
else
   touch results.tsv
fi
errorcounter=0
echo -e "Number of Isoforms \t Found isoforms by IsONform \n">>results.tsv
ls
file=out/mapping.txt
echo "this is the whole list of dir"
for ((i=3; i<=50; i++))
do
n_reads=$(($i+5))
#echo "Generating $i TestIsoforms" >>results.tsv
for((j=1;j<=5;j++))
do
outputs=0
python generateTestCases.py --ref /home/alexanderpetri/Desktop/RAWDATA_PhD1/Isoform_Test_data.fa --sim_genome_len 1344 --nr_reads $n_reads --outfolder testout --coords 50 100 150 200 250 300 350 400 450 500 --probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms $i --e False
python main.py --fastq ~/PHDProject1/testout/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout 
res=$(< "$file" wc -l)
result=$(($res/2))
echo "$result"
echo -e "$i \t $result \t" >> results.tsv
if [[ "$i" != "$result" ]]
then
echo "Error found!"
cp testout/isoforms.fa error_${errorcounter}.fa
errorcounter=$((errorcounter+1))
echo "New file created" >>results.tsv
fi
done
echo -e"\n">>results.tsv
done

