QMC=Chem : Quantum Monte Carlo for Chemistry
============================================

QMC=Chem is the quantum Monte Carlo program of the
[Toulouse (France) group](http://qmcchem.ups-tlse.fr).
It is meant to be used in the *post-Full-CI* context : a quasi-Full-CI
calculation is done with the
[quantum package](https://github.com/LCPQ/quantum_package),
and this wave function is used as a trial wave function for the fixed-node
diffusion Monte Carlo algorithm.


* Parallel efficiency of 98.3% on 16_000 cores
* The load balancing is optimal: the workers always work 100% of the time,
  independently of their respective CPU speeds
* Efficient: 0.96 Pflops/s on 76_800 cores of Curie in 2011
* All network communications are non-blocking,
  with the [ZeroMQ](http://zeromq.org) library
* All the implemented algorithms are CPU-bound : the only limit
  is the available CPU time
* The number of simultaneous worker nodes can be variable during a calculation
* Fully fault-tolerant (crashing nodes don't stop the running calculation)
* QMC=Chem has been used in grid environments (EGI european grid) and 
  in Cloud environments (France Grilles) coupled to supercomputers 


Warnings:
* QMC=Chem is under the GPLv2 license. Any modifications to or
  software including (via compiler) GPL-licensed code must also be made available
  under the GPL along with build & install instructions.


Requirements
------------

* [OCaml compiler with Opam and Core library](http://github.com/ocaml)
* [ZeroMQ high performance communication library](http://www.zeromq.org)
* [F77_ZMQ ZeroMQ Fortran interface](http://github.com/scemama/f77_zmq/)
* [IRPF90 Fortran code generator](http://irpf90.ups-tlse.fr)
* [EZFIO Easy Fortran I/O library generator](http://github.com/scemama/EZFIO)
* GNU C++ Compiler (g++) for ZeroMQ 
* Python >= 2.6 for install scripts
* Bash for install scripts
* Fortran compiler, Intel Fortran recommended
* Lapack library, Intel MKL recommended


Most of the dependencies are open-source can be downloaded automatically by
going into the `install` directory and running `make`. It will first download
into the `install/Downloads` directory everything that needs to be installed.
The building of the dependencies takes place in the `install/_build`
directory, and the packages that are being installed can be followed by looking
at the log files in this directory. When a package was successfully installed,
a `*.ok` file is created and the log file is deleted.

If you don't have an internet connection available, you can execute the
downloading step on another computer and transfer all the downloaded files
into the `Downloads` directory.
The Fortran and C++ compilers, Python and Bash interpreters and the Lapack
library need to be installed manually by the user.

Installation
------------

The `make.config` file contains compiler specific parameters. You should change
them to match your hardware. You can copy the `make.config.ifort` or
`make.config.gfortran` as a starting point.

Before using or compiling QMC=Chem, environment variables need to be loaded. The
environment variables are located in the `qmcchemrc` file:

```bash
$ source qmcchemrc
```

The `QMCCHEM_NIC` environment variable should be set to the proper network interface,
usually `ib0` on HPC machines.

To compile the program, run

```bash
$ make
```


Example of a QMC=Chem calculation
---------------------------------

1.Calculation with the [quantum package](http://github.com/QuantumPackage/qp2)


Create the `xyz` file containing the nuclear coordinates of the system

```bash
$ cat > h2o.xyz << EOF
3
Water molecule
O 0. 0. 0.
H 0.9572 0. 0.
H -0.239987 0.926627 0.
EOF
```

Choose a suitable basis set and create the [EZFIO database](https://github.com/LCPQ/ezfio)

```bash
$ qp_create_ezfio -b cc-pvdz h2o.xyz -o h2o
```

Run the SCF calculation

```bash
$ qp_run scf h2o
```
Run the CIPSI calculation

```bash
$ qp_run fci h2o
```

Transform the input for use in QMC=Chem

```bash
$ qp_run save_for_qmcchem h2o
```


2.FN-DMC calculation with QMC=Chem


Before using QMC=Chem, you need to load the environment variables:

```bash
$ source qmcchemrc
```

In QMC=Chem, everything goes through the use of the ``qmcchem`` command.
When a command is run with no arguments, it prints a help message.
This is mainly the manual of QMC=Chem. For example:

```
$ qmcchem 
QMC=Chem command

  qmcchem SUBCOMMAND

=== subcommands ===

  debug    Debug ZeroMQ communications
  edit     Edit input data
  md5      Manipulate input MD5 keys
  result   Displays the results computed in an EZFIO directory.
  run      Run a calculation
  stop     Stop a running calculation
  version  print version information
  help     explain a given subcommand (perhaps recursively)

missing subcommand for command qmcchem

$ qmcchem edit
Run a calculation

  qmcchem run EZFIO_FILE


Run QMC=Chem
      

=== flags ===

  [-a]                    Add more resources to a running calculation.
  [-d]                    Start a dataserver process on the local host.
  [-q <dataserver_addr>]  Start a qmc process on the local host.
  [-s <host>]             Start a qmc process on <host>.
  [-help]                 print this help text and exit
                          (alias: -?)

missing anonymous argument: EZFIO_FILE
```

1) Set the parameters for a VMC calculation to create initial walker positions

```
$ qmcchem edit -h
Edit input data

  qmcchem edit EZFIO_FILE [INPUT]


Edit input data
      

=== flags ===

  [-c]                Clear blocks
  [-e energy]         Fixed reference energy to normalize DMC weights
  [-f 0|1]            Correct wave function to verify electron-nucleus cusp
                      condition
  [-j jastrow_type]   Type of Jastrow factor [ None | Core | Simple ]
  [-l seconds]        Length (seconds) of a block
  [-m method]         QMC Method : [ VMC | DMC ]
  [-n norm]           Truncation t of the wave function : Remove determinants
                      with a
                      contribution to the norm less than t
  [-s sampling]       Sampling algorithm : [ Langevin | Brownian ]
  [-t seconds]        Requested simulation time (seconds)
  [-ts time_step]     Simulation time step
  [-w walk_num]       Number of walkers per CPU core
  [-wt walk_num_tot]  Total number of stored walkers for restart
  [-help]             print this help text and exit
                      (alias: -?)

$ qmcchem edit h2o -f 1 -m VMC -n 1.e-5 -s Langevin -t 300 -l 10
```

3) Get info on the wave function

```bash
$ qmcchem info h2o
```

4) Run the VMC calculation

```bash
$ qmcchem run h2o
```

5) Set the correct parameters for FN-DMC

```bash
$ qmcchem edit h2o -e -76.438 -m DMC -s Brownian -ts 3.e-4 -t 3600 -l 30
```

6) Run the FN-DMC calculation

```bash
$ qmcchem run h2o
```

7) Print the result

```bash
$ qmcchem result h2o

```




References
----------

[Quantum Monte Carlo for large chemical systems: Implementing efficient strategies for petascale platforms and beyond](http://dx.doi.org/10.1002/jcc.23216)
> Anthony Scemama , Michel Caffarel , Emmanuel Oseret and William Jalby (2013), in: Journal of Computational Chemistry, 34:11(938--951) 

[Quantum Monte Carlo with very large multideterminant wavefunctions](http://arxiv.org/abs/1510.00730)
> Anthony Scemama , Thomas Applencourt , Emmanuel Giner and Michel Caffarel (2015), in: ArXiv ePrints:arXiv:1510.00730v2 [physics.chem-ph] 

