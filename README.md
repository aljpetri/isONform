# IsONform(er)
The first project for my PHD:
This tool generates isoforms out of clustered and corrected ONT reads.
For this a graph is built up using the networkx api and different simplification strategies are applied to it.
The algorithm uses spoa to generate the final isoforms.It produces two files:
-mapping.txt contains information about which reads were mapped together into which consensus. It has the following form
Line1:consensusID
Line2: List of read names

-spoa.fa contains the actual isoforms stored in the fasta format:
Line1: >consensusID
Line2: consensus sequence

To run the code:
python main.py --fastq /home/alexanderpetri/Desktop/RAWDATA_PhD1/RBMYnoDeletionsNoSC.fa --k 9 --w 10 --xmin 14 --xmax 80  --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder out

