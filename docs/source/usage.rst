.. _basic-usage-and-examples-qmc-chem:

Basic Usage and Examples of QMC=Chem
====================================

This section provides a comprehensive guide on how to perform basic operations and run simple examples using QMC=Chem.

Example of a QMC=Chem Calculation
---------------------------------

Preparing the System with Quantum Package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Create the EZFIO Database**

   Begin by creating an ``xyz`` file with the nuclear coordinates of the system. For a water molecule:

   .. code-block:: bash

       $ cat > h2o.xyz << EOF
       3
       Water molecule
       O 0. 0. 0.
       H 0.9572 0. 0.
       H -0.239987 0.926627 0.
       EOF

2. **Select Basis Set and Pseudopotential**

   Choose an appropriate basis set and pseudopotential, then create the EZFIO database:

   .. code-block:: bash

       $ qp_create_ezfio -b cc-pvdz_ecp_bfd -p bfd h2o.xyz -o h2o

3. **Run SCF and CIPSI Calculations**

   Perform the SCF calculation, followed by a CIPSI calculation:

   .. code-block:: bash

       $ qp_run scf h2o
       $ qp_run fci h2o

4. **Transform Input for QMC=Chem**

   Prepare the input for QMC=Chem:

   .. code-block:: bash

       $ qp_run save_for_qmcchem h2o

Performing FN-DMC Calculation with QMC=Chem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Load Environment Variables**

   Before starting, load the QMC=Chem environment:

   .. code-block:: bash

       $ source qmcchemrc

2. **Understanding the ``qmcchem`` Command**

   The ``qmcchem`` command is central to operations in QMC=Chem. Running it without arguments displays a help message:

   .. code-block:: none

       $ qmcchem
       [Output showing various subcommands]

3. **Set Parameters for VMC Calculation**

   Configure the VMC calculation to create initial walker positions:

   .. code-block:: none

       $ qmcchem edit h2o -f 1 -m VMC -n 1.e-5 -s Langevin -t 300 -l 10

4. **Run the VMC Calculation**

   Execute the VMC calculation:

   .. code-block:: bash

       $ qmcchem run h2o

5. **Configure FN-DMC Parameters**

   Adjust the parameters for the FN-DMC calculation:

   .. code-block:: bash

       $ qmcchem edit h2o -e -76.438 -m DMC -s Brownian -ts 3.e-4 -t 3600 -l 30

6. **Execute FN-DMC Calculation**

   Perform the FN-DMC calculation:

   .. code-block:: bash

       $ qmcchem run h2o

7. **View Results**

   Finally, view the results of the calculation:

   .. code-block:: bash

       $ qmcchem result h2o

