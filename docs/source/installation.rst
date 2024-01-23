.. _installation-and-configuration:

Installation and Configuration
==============================

This section provides detailed instructions for installing and configuring
QMC=Chem, including prerequisites, supported platforms, and initial setup
steps.

Requirements
------------

Before installing QMC=Chem, ensure that the following prerequisites are met:

- **Bash**: A Unix shell and command language.
- **Python3**: Required for various scripting and automation tasks.
- **Fortran Compiler**: Intel Fortran is recommended for optimal performance.
- **Lapack Library**: Intel MKL is recommended for efficient numerical computations.
- **ZeroMQ**: A high-performance communication library, essential for network communications. `ZeroMQ Library <http://www.zeromq.org>`_
- **F77_ZMQ**: A ZeroMQ Fortran interface. `Download F77_ZMQ <https://github.com/zeromq/f77_zmq/releases/download/v4.3.3/f77-zmq-4.3.3.tar.gz>`_
- **QMCkl Library**: A library for quantum Monte Carlo kernels. `Download QMCkl <https://github.com/TREX-CoE/qmckl/releases/download/v0.3.1/qmckl-0.3.1.tar.gz>`_
- **TREXIO Library**: A library for managing quantum chemistry data. `Download TREXIO <https://github.com/TREX-CoE/trexio/releases/download/v2.2.0/trexio-2.2.0.tar.gz>`_
- **OCaml Compiler with Opam**: Required for certain computational tasks. `OCaml Compiler <http://github.com/ocaml>`_

To install the necessary OCaml packages, execute the following command:

.. code-block:: bash

    opam install ocamlbuild cryptokit zmq sexplib ppx_sexp_conv ppx_deriving getopt trexio

For OCaml installation issues on x86 systems, consider downloading
`this archive <https://github.com/QuantumPackage/qp2-dependencies/blob/e0d0e02e9f5ece138d1520106954a881ab0b8db2/x86_64/opampack.tar.gz>`_
and follow these steps:

.. code-block:: bash

    tar --gunzip --extract --file opampack.tar.gz
    cd opampack
    ./install.sh
    export OPAMROOT="${PWD}"/opamroot
    eval ("${OPAMROOT}"/opam env)

Installation Process
--------------------

To compile QMC=Chem, execute the following commands:

.. code-block:: bash

    ./autogen.sh
    ./configure && make

After compilation, it's necessary to set up the environment for QMC=Chem. The environment variables are specified in the ``qmcchemrc`` file:

.. code-block:: bash

    source qmcchemrc

Set the ``QMCCHEM_NIC`` environment variable to the appropriate network interface, typically ``ib0`` on HPC systems.
If needed, set ``QMCCHEM_IO=b`` to store the results in binary format to accelerate post-processing.

Tuning for HPC Systems
----------------------

For optimal performance, the Intel compiler is recommended. If using gfortran version 12 or higher, add the ``-fallow-argument-mismatch`` option.

Example configurations:

- **For Intel Fortran Compiler**:

  .. code-block:: none

      ./configure FC=ifort FCFLAGS="-O2 -g -ip -ftz -finline -xCORE-AVX2 -mkl=sequential"

- **For GCC version 12 or higher**:

  .. code-block:: none

      ./configure FCFLAGS="-fallow-argument-mismatch -g -O2 -ffast-math -march=native -fno-trapping-math -fno-math-errno -ftree-vectorize -fno-stack-protector -fopenmp"

File Preparation
----------------

To prepare files for QMC=Chem, the ``save_for_qmcchem`` plugin must be installed in Quantum Package:

.. code-block:: bash

    qp plugins download https://gitlab.com/scemama/qp_plugins_scemama
    qp plugins install qmcchem
    cd $QP_ROOT/src/qmcchem
    ninja

After completing a Quantum Package calculation, execute the following command to prepare the directory for QMC=Chem:

.. code-block:: bash

    qp run save_for_qmcchem

