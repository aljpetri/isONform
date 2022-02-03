# isONform- an algorithm capable of recovering isoforms from long read sequencing data
# Table of contents
1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Output](#output) 
4. [Running the test script](#Running)
	1. [Running isONform](#runalgo)
5. [Credits](#credits)

## Installation <a name="installation"></a>

TBD

### Dependencies

1. `networkx`
2. `ordered-set`
3. `matplotlib`
4. `parasail`
5. `edlib`


## Introduction <a name="introduction"></a>

This tool generates isoforms out of clustered and corrected long reads.
For this a graph is built up using the networkx api and different simplification strategies are applied to it, such as bubble popping and node merging.
The algorithm uses spoa to generate the final isoforms.<br />

## Output <a name="output"></a>

The algorithm produces two files:<br />
-<strong>mapping.txt</strong> contains information about which reads were mapped together into which consensus. It has the following form:<br />
Line1:consensusID <br />
Line2: List of read names </p>

-<strong>spoa.fa</strong> contains the actual isoforms stored in the fasta format:<br />
Line1: >consensusID<br />
Line2: consensus sequence<br />
## Running the code <a name="Running"></a>
If you want to generate Simulated Isoforms for testing:<br />
(On my machine:)<br />
```
python generateTestCases.py --ref /home/alexanderpetri/Desktop/RAWDATA_PhD1/Isoform_Test_data.fa 
					--sim_genome_len 1344 --nr_reads 10 --outfolder testout 
					--coords 50 100 150 200 250 300 350 400 450 500 
					--probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms 8 
```

```
python generateTestCases.py --ref /path/to/Isoform_Test_data.fa 
							--sim_genome_len 1344 --nr_reads 10 --outfolder testout 
							--coords 50 100 150 200 250 300 350 400 450 500 
							--probs 0.4 0.4 0.4 0.4 0.4 --n_isoforms 8
```

### Actual algorithm <a name="runalgo"></a>
To run the actual algorithm:<br />
(On my machine:)

```
python main.py --fastq ~/PHDProject1/testout/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout
```

```
python main.py --fastq /path/to/isoforms.fa --k 9 --w 10 --xmin 14 --xmax 80 --exact --max_seqs_to_spoa 200 --max_bubblesize 2 --delta_len 3 --outfolder testout
```

## Credits <a name="credits"></a>

Please cite [1] when using isONform.
