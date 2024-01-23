![QMC=Chem logo](https://github.com/trex-coe/qmcchem2/raw/master/docs/source/_static/QmcChemLogo.png)

This repository contains version 2 of QMC=Chem.
Version 1 is available on [GitLab](https://gitlab.com/scemama/qmcchem).
QMC=Chem is a quantum Monte Carlo program meant to be used after
preparing a trial wave function with the
[Quantum Package](https://github.com/quantumpackage/qp2) code.

The documentation is available on [ReadTheDocs](https://qmcchem2.readthedocs.io/en/latest/index.html).

Requirements
------------

* Bash
* Python3
* Fortran compiler, Intel Fortran recommended
* Lapack library, Intel MKL recommended
* [ZeroMQ high performance communication library](http://www.zeromq.org)
* [F77_ZMQ ZeroMQ Fortran interface](http://github.com/zeromq/f77_zmq/): https://github.com/zeromq/f77_zmq/releases/download/v4.3.3/f77-zmq-4.3.3.tar.gz
* [QMCkl library](https://github.com/trex-coe/qmckl): https://github.com/TREX-CoE/qmckl/releases/download/v0.3.1/qmckl-0.3.1.tar.gz
* [TREXIO library](https://github.com/trex-coe/trexio): https://github.com/TREX-CoE/trexio/releases/download/v2.2.0/trexio-2.2.0.tar.gz
* [OCaml compiler with Opam](http://github.com/ocaml)

To install the required OCaml packages, run
```bash
opam install ocamlbuild cryptokit zmq sexplib ppx_sexp_conv ppx_deriving getopt trexio
```

If you have trouble installing OCaml, on x86 systems you can download
the [this file](https://github.com/QuantumPackage/qp2-dependencies/raw/master/ocaml-bundle_x86.tar.gz)
and run
```bash
tar -zxf ocaml-bundle_x86.tar.gz
./ocaml-bundle/bootstrap.sh
./ocaml-bundle/configure.sh
./ocaml-bundle/compile.sh 
```


Installation
------------

To compile the program, run

```bash
$ ./autogen.sh
$ ./configure && make
```

Before using QMC=Chem, environment variables need to be loaded. The
environment variables are located in the `qmcchemrc` file:

```bash
$ source qmcchemrc
```

The `QMCCHEM_NIC` environment variable should be set to the proper network interface,
usually `ib0` on HPC machines.

Tuning for the HPC system
-------------------------

The Intel compiler is preferred for better performance.
With gfortan >= 12, the option `-fallow-argument-mismatch` should be added.

Here are a few examples:
```
# Ifort
./configure FC=ifort FCFLAGS="-O2 -g -ip -ftz -finline -xCORE-AVX2 -mkl=sequential"

# GCC >= 12
./configure FCFLAGS="-fallow-argument-mismatch -g -O2  -ffast-math  -march=native   -fno-trapping-math -fno-math-errno -ftree-vectorize -fno-stack-protector -fopenmp"
```


File preparation
----------------

To create files suitable for QMC=Chem, the `save_for_qmcchem` plugin
needs to be installed in Quandtum Package. This can be done as
```bash
qp plugins download https://gitlab.com/scemama/qp_plugins_scemama
qp plugins install qmcchem
cd $QP_ROOT/src/qmcchem
ninja
```

Then, after running a Quantum Package calculation you should run
```bash
qp run save_for_qmcchem
```
to prepare the directory for use in QMC=Chem.




-----------------
![European flag](https://trex-coe.eu/sites/default/files/inline-images/euflag.jpg)
[TREX: Targeting Real Chemical Accuracy at the Exascale](https://trex-coe.eu) project has received funding from the European Union’s Horizon 2020 - Research and Innovation program - under grant agreement no. 952165. The content of this document does not represent the opinion of the European Union, and the European Union is not responsible for any use that might be made of such content.
