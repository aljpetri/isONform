# isONform - Reference-free isoform reconstruction from long read sequencing data
# Table of contents
1. [Installation](#installation)
2. [Introduction](#introduction)
3. [Output](#output) 
4. [Input data](#Input_data)
5. [Running isONform](#Running)
	1. [Running a test](#runtest)
6. [Credits](#credits)

## Installation <a name="installation"></a>

isONform has been re-implemented in Rust (2026-09-03). It needs a Rust toolchain
([rustup.rs](https://rustup.rs)).

```
git clone https://github.com/aljpetri/isONform.git
cd isONform/rust
cargo build --release
```

That produces `target/release/isONform_parallel` and `target/release/main`. Put
them on your `PATH` and run them as shown under
[Running isONform](#Running).

The original python implementation is still available and is the reference this
port is checked against; see [INSTALL-python.md](INSTALL-python.md).

### Which version of the rust-port

By default, the rust port uses the WFA2 aligner, and is much faster and
slightly more accurate on average, but does not produce identical results to the Python version. 
The **`--faithful`** parameter reproduces the python implementation **byte for byte**, but
is only about **~2x faster than python** on shallow data and about the
same speed on the deepest clusters. We recommend using the port in default mode (no `--faithful` flag). 


Full comparison against the python implementation --- accuracy,
redundancy, runtime and peak memory on five corpora from 10 000 to 1 000 000 reads
are found here: [Port-benchmark.md](Port-benchmark.md).

### Running a test <a name="runtest"></a>

`test_data/sirv_sim` is a small dataset for checking an installation: two
clusters of 51 corrected simulated SIRV reads each.

```
isONform_parallel --fastq_folder test_data/sirv_sim --outfolder /tmp/isonform_test \
                  --t 4 --split_wrt_batches --iso_abundance 3
```

This should finish in seconds and write **7 isoforms** to
`/tmp/isonform_test/transcriptome.fasta`, recovering SIRV101, SIRV102, SIRV103,
SIRV201, SIRV202 and SIRV204 at >= 0.993 identity. See
[test_data/README.md](test_data/README.md) for provenance and for what the
numbers mean.

## Introduction <a name="introduction"></a>

IsONform generates isoforms out of clustered and corrected long reads.
For this a graph is built up using the networkx api and different simplification strategies are applied to it, such as bubble popping and node merging.
The algorithm uses spoa to generate the final isoforms.<br />
## Input data <a name="Input_data"></a>
The isONpipeline takes .fastq files generated with long-read sequencing techniques (ONT or Pacbio) as an input that additionally have been cleaned of barcodes.
Please make sure that you run the isONpipeline on data that have been processed with  [LIMA](https://lima.how/) (Pacbio data) or [Pychopper](https://github.com/epi2me-labs/pychopper) (ONT data) so that all the barcodes are removed from the reads

## Running isONform <a name="Running"></a>

To only run the isONform algorithm:<br />


```
isONform_parallel --fastq_folder path/to/input/files --t <nr_cores> --outfolder /path/to/outfolder --split_wrt_batches 
```

Argument names, defaults, validation messages and exit codes match the python
implementation, so any existing command or script works unchanged. Add
`--faithful` to reproduce the python output byte for byte.

The full isON-pipeline (isONclust, isONcorrect, isONform) can be found [here](https://github.com/aljpetri/isONform/blob/master/isON_pipeline.sh) and is run via:

```
./isON_pipeline.sh --raw_reads </absolute/path/to/raw_reads.fq>  --outfolder <outfolder>  --num_cores <num_cores> --isONform_folder <isONform_folder> --iso_abundance <iso_abundance> --mode <mode>
```
(Please note that this requires isONclust [LINK](https://github.com/ksahlin/isONclust) and isONcorrect [LINK](https://github.com/ksahlin/isONcorrect) to be installed in addition to isONform)

To receive more information about the arguments used for the isON_pipeline script:
```
./isON_pipeline.sh --help
```

## Outputs <a name="Outputs"></a>
IsONform outputs three main files: transcriptome.fasta, mapping.txt, and support.txt.
For each isoform that isONform reconstructs the id has the following form: x_y_z.

'x' denotes the isONclust cluster that the isoform stems from.
As we cluster reads as in isONcorrect in batches of 1000 reads the 'y' denotes from which batch the isoform was reconstructed.
The 'z' denotes a unique identifier which enables us to have unique ids for each isoform that we reconstructed.
In mapping.txt it is indicated from which original reads an isoform has been reconstructed.
support_txt gives the support (i.e. how many original reads make up the isoform).

## Contact <a name="Contact"></a>
If you encounter any problems, please raise an issue on the issues page.

## Credits <a name="credits"></a>

Please cite [1] when using isONform.

1. Petri, A. J., & Sahlin, K. (2023). isONform: reference-free transcriptome reconstruction from Oxford Nanopore data. Bioinformatics, 39(Supplement_1), i222-i231. https://academic.oup.com/bioinformatics/article/39/Supplement_1/i222/7210488 .

Please additionally cite [2] and [3] when running the full pipeline.

2. Kristoffer Sahlin, Paul Medvedev. De Novo Clustering of Long-Read Transcriptome Data Using a Greedy, Quality-Value Based Algorithm, Journal of Computational Biology 2020, 27:4, 472-484. [Link](https://www.liebertpub.com/doi/abs/10.1089/cmb.2019.0299).
3. Sahlin, K., Medvedev, P. Error correction enables use of Oxford Nanopore technology for reference-free transcriptome analysis. Nat Commun 12, 2 (2021). https://doi.org/10.1038/s41467-020-20340-8  [Link](https://www.nature.com/articles/s41467-020-20340-8).
