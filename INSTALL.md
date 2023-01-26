# Installation instructions for QMC=Chem

1. Install TREXIO and QMCKL external libraries.

2. Download the QMC=Chem source package.

3. Extract the package to a folder:
```
$ tar -xzf qmcchem-2.0.0.tar.gz
```

4. Change to the directory of the extracted source package:
```
$ cd qmcchem-2.0.0
```

5. Configure the build with:
```
$ ./configure
```
The default compiler is the Intel Fortran compiler (`ifort`) if detected,
otherwise GNU Fortran (`gfortran`) will be used.
If the Intel compiler is detected, the MKL library is used, otherwise the
user's BLAS library is checked.

6. Compile the source code with:
```
$ make
```

7. Install the program with:
```
$ make install
```

Optional configuration arguments:

* `--with-arch`: Specifies the target architecture to generate code for. Can be `avx`, `avx2` or `avx512`. The default is `avx2`.

-----

To install QMC=Chem with the `--with-arch` option, use the following commands:

```
./configure --with-arch=<avx|avx2|avx512> && make
```



