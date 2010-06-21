
.. _mira-assembly:

Mira assembly
-------------

The `mira <http://sourceforge.net/apps/mediawiki/mira-assembler/index.php?title=Main_Page>`_ assembler is used to create a set of contigs with the sequence reads. The reads can come from 454, sanger and illumina sequencing. Hybrid assembly are possible. For mira configuration details refer to its documentation.

Input and output files
______________________

The input files required to do a mira analysis are the reads located in reads/cleaned. The reads files should follow the :doc:`naming conventions <introduction>`.

Configuration parameters
________________________

The default configuration is tailored to EST assemblies. To modify the mira command line parameters you should go to the mira section in the ngs_backbone.conf file. The options are:

job_options
  The mira job options parameter. By default they are: denovo, est

general settings
  The parameters that affect all platforms.

454_settings
  The parameters that affect the 454 reads.

sanger_settings
  The parameters that affect the sanger reads.

Running the analysis
____________________

The mira assembly analysis is divided into three ngs_backbone analyses: prepare_mira_assembly, mira_assembly and select_last_assembly.

The analysis prepare_mira_assembly will create the files required as input by mira in the directory assembly/input/. These files will be created taking the reads from reads/cleaned/.

The mira_assembly analysis runs mira and creates the contigs. The files created by this analysis will be located at a timestamped directory located in assembly/. Several assemblies could be created with different parameters and each one would go into a different timestamped directory. Inside these directories a result subdirectory is created with the relevant result files.

The select_last_assembly will just make a soft link named assembly/result that points to the result subdirectory located in the latest timestamped assembly.


.. _mapping:
