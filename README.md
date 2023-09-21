# isONform - Reference-free isoform reconstruction from long read sequencing data
# Table of contents
1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Output](#output) 
4. [Running isONform](#Running)
	1. [Running a test](#runtest)
5. [Credits](#credits)

## Installation <a name="installation"></a>


### Via pip
```
pip install isONform
```

This command installs isONforms dependencies:

1. `networkx`
2. `ordered-set`
3. `matplotlib`
4. `parasail`
5. `edlib`
6. `pyinstrument`
7. `namedtuple`
8. `recordclass`


### From github source
1. Create a new environment for isONform (at least python 3.7 required):<br />
		`conda create -n isonform python=3.10 pip` <br />
		`conda activate isonform` <br />
2.  Install isONcorrect and SPOA <br />
		`pip install isONcorrect` <br />
		`conda install -c bioconda spoa` <br />
3.  Install other dependencies of isONform:<br />
		`conda install networkx`<br />
		`pip install parasail`<br />

4. clone this repository


## Introduction <a name="introduction"></a>

IsONform generates isoforms out of clustered and corrected long reads.
For this a graph is built up using the networkx api and different simplification strategies are applied to it, such as bubble popping and node merging.
The algorithm uses spoa to generate the final isoforms.<br />


## Running isONform <a name="Running"></a>

To run the algorithm:<br />


```
python isONform_parallel.py --fastq_folder path/to/input/files --t <nr_cores> --outfolder /path/to/outfolder --split_wrt_batches 
```

the isON-pipeline (isONclust, isONcorrect, isONform) can be run via:

```
./full_pipeline.sh <raw_reads.fq>  <outfolder>  <num_cores> <isONform_folder> <iso_abundance> <mode>
```
(Note that this requires pychopper, isONclust and isONcorrect to be installed)


## Contact <a name="Contact"></a>
If you encounter any problems, please raise an issue on the issues page, you can also contact the developer of this repository via:
alexander.petri[at]math.su.se


## Credits <a name="credits"></a>

Please cite [1] when using isONform.

1. Petri, A. J., & Sahlin, K. (2023). isONform: reference-free transcriptome reconstruction from Oxford Nanopore data. Bioinformatics, 39(Supplement_1), i222-i231. https://academic.oup.com/bioinformatics/article/39/Supplement_1/i222/7210488 .

