
.. _mira-assembly:

Mira assembly
-------------

The mira_ assembler can be used to create a set of contigs with the sequence reads.
The reads can come from 454, sanger and illumina sequencing.
Hybrid assemblies are possible.
For mira configuration details refer to its documentation.

Running the analysis
____________________

The input files required to do a mira analysis are the reads located in reads/cleaned. The reads files should follow the :doc:`naming conventions <introduction>`.

The mira assembly analysis prepares the files required by mira.

The analysis prepare_mira_assembly will create the files required as input by mira in the directory assembly/input/. These files will be created taking the reads from reads/cleaned/.

.. include:: links.txt

