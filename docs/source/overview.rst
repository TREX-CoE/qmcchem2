.. _overview:

Overview
========

Technical Features and Computing Performance
--------------------------------------------

The program exhibits good parallel efficiency, achieving approximately 98.3% on
16,000 cores. This efficiency reflects the program's ability to utilize
computational resources effectively. Load balancing is a key feature of
QMC=Chem, ensuring consistent workload distribution across different CPU
speeds.

In terms of computational achievements, QMC=Chem reached 0.96 Pflops/s on
76,800 cores of the Curie supercomputer in 2011. This performance indicates the
program's capability to scale up in high-performance computing environments.
The use of non-blocking network communications via the ZeroMQ library aids in
maintaining efficient data transfer.

QMC=Chem's algorithms are primarily CPU-bound, meaning the program's
performance scales with the available CPU time. This aspect allows for
flexibility in computational tasks, depending on the available resources. The
program also supports a variable number of worker nodes during calculations,
which can be beneficial for managing computational loads dynamically.


Reliability and Adaptability
----------------------------

The program is designed with fault tolerance in mind, allowing calculations to
continue despite potential node failures. This feature is particularly useful
in less predictable computing environments like grids or clouds.

QMC=Chem has been utilized in various settings, including grid environments
like the EGI European grid and cloud platforms such as France Grilles, in
combination with supercomputers. This demonstrates its adaptability to
different computational infrastructures.


Licensing and Community Engagement
----------------------------------

QMC=Chem is available under the GPLv2 license, encouraging open-source
collaboration. This licensing means that any modifications or derivative works
that include GPL-licensed code are also required to be open-source, along with
the necessary build and installation instructions. This approach supports a
collaborative development environment, allowing for shared improvements and
enhancements to the program.


