# Example of a QMC=Chem calculation

## Calculation with [quantum package](http://github.com/QuantumPackage/qp2)


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

Choose a suitable basis set and pseudopotential and create the [EZFIO database](https://gitlab.com/scemama/ezfio)

```bash
$ qp_create_ezfio -b cc-pvdz_ecp_bfd -p bfd h2o.xyz -o h2o
```

Run the SCF calculation

```bash
$ qp_run scf h2o
```

Run a CIPSI calculation

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

