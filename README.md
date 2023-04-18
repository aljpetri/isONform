# isONform- an algorithm capable of recovering isoforms from long read sequencing data
# Table of contents
1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Output](#output) 
4. [Running isONform](#Running)
	1. [Running a test](#runtest)
5. [Credits](#credits)

## Installation <a name="installation"></a>

### Dependencies

1. `networkx`
2. `ordered-set`
3. `matplotlib`
4. `parasail`
5. `edlib`
6. `pyinstrument`
7. `namedtuple`
8. `recordclass`


## Installation guide
### Via pip
```
pip install isONform
```

### Downloading source from github
1. Create a new environment for isONform (at least python 3.7 required):<br />
		`conda create -n isonform python=3.10 pip` <br />
		`conda activate isonform` <br />
2.  Install isONcorrect and SPOA (strongly recommended)<br />
		`pip install isONcorrect` <br />
		`conda install -c bioconda spoa` <br />
3.  Install other dependencies of isONform:<br />
		`conda install networkx`<br />
		`pip install ordered-set`<br />
		`conda install matplotlib`<br />
		`pip install parasail`<br />
		`pip install pyinstrument`<br />
		`conda install -c cerebis recordclass`<br />
4. clone this repository
		




## Introduction <a name="introduction"></a>

This tool generates isoforms out of clustered and corrected long reads.
For this a graph is built up using the networkx api and different simplification strategies are applied to it, such as bubble popping and node merging.
The algorithm uses spoa to generate the final isoforms.<br />

## Output <a name="output"></a>

The algorithm produces three files:<br />
-<strong>mapping.txt</strong> contains information about which reads were mapped together into which consensus. It has the following form:<br />
Line1:consensusID <br />
Line2: List of read names </p>

-<strong>transcriptome.fastq</strong> contains the actual isoforms stored in the fastq format:<br />

-<strong>support.txt</strong> contains the consensusID as well as the number of reads supporting it<br/>
## Running isONform <a name="Running"></a>

To run the algorithm:<br />


```
python isONform_parallel.py --fastq_folder path/to/input/files --t <nr_cores> --exact_instance_limit 50 --k 20 --w 31 --xmin 14 --xmax 80 --max_seqs_to_spoa 200 --delta_len 5 --outfolder /path/to/outfolder --iso_abundance 3 --split_wrt_batches --delta_iso_len_3 5 --delta_iso_len_5 5
```

the isON-pipeline can be run via:
```
./full_pipeline.sh    <raw_reads.fq>  <outfolder>  <num_cores> <isONform_folder> <iso_abundance> <mode>
```
(Note that this requires pychopper,isONclust and isONcorrect to be installed)



### Running the test analysis pipeline <a name="runtest"></a>

To run the test analysis pipeline:

```
./generateTestResults.sh  </path/to/input/reference.fa> <output_root>
```

If you want to generate Simulated Isoforms for testing

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


## Credits <a name="credits"></a>

Please cite [1] when using isONform.
