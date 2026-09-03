# Installing the python implementation of isONform

The Rust implementation is the recommended one and needs only a Rust toolchain;
see [Installation](README.md#installation) in the README.

These are the instructions for the original python implementation, which is still
the reference the Rust one is verified against.

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
