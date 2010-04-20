
Usage
=====

The main ideas to consider when using franklin are: project and analysis. A *project* is a directory (with its subdirectories) that includes a configuration file and all input and output files. An *analysis* takes some inputs from the project and creates some outputs. franklin knows where the input and output files are for every analysis because the project directory structure is the same for every project. For instance the reads to clean are always in the directory /reads/original and the cleaned reads always are in /reads/cleaned/.

The configuration parameters required for every analysis are stored in the configuration file (franklin.conf) located in the project directory. A sample of this configuration file is created by franklin when a new project is created. This file should be tweaked to adapt the analyses to your requirements before running them.

franklin has only two main executables: franklin_create_project.py and franklin_analysis.py. The first one is used to create a project from scratch and the second one to run all analyses.

Creating a new project
----------------------

::

  $ franklin_create_project.py -p project_name
  $ ls -l project_name
  total 4
  -rw-r--r-- 1 jose jose 2479 abr 16 09:30 franklin.conf

This command will create a new directory named project_name with a file named franklin.conf in it. These are the two hallmarks that define a franklin *project*, the directory and the configuration file.

Running an analysis
-------------------

::

  $ franklin_analysis.py -a analysis_name

This command will run the analysis on the data present in the project using the parameters found in the configuration file. The output files will also be located in the project.

Available analyses
------------------

The available analyses are:

========================================    =================================================
analysis                                    description
========================================    =================================================
:ref:`clean-reads`                          sequence reads cleaning
:ref:`mira-assembly`                        Assembly reads into a contig set with  mira
mapping                                     bwa read mapping against a reference sequence
bam realignment                             GATK bam realignment
SNP calling                                 SNP annotation from a bam file
ORF annotation                              ESTScan ORF annotation
ortholog annotation                         reciprocal blast based ortholog annotation
description annotation                      description blast based annotation
SSR annotation                              microsatellite sputnik based annotation
cdna intron annotation                      est2genome cDNA based annotation
GO annotation                               blast2go  annotation
========================================    =================================================

Naming conventions
------------------

The franklin usage is heavily based on directory and file name conventions. If the input files are not located where franklin expects to find them or they have non-standard file names the analysis will fail. 

====================== ========================
file or directory type location
====================== ========================
configuration file     franklin.conf
log file               franklin.log
raw reads              /reads/original/
clean reads            /reads/cleaned/
assemblies             /assembly/
assembly input         /assembly/input/
assembly output        /assembly/result/
mappings               /mapping/
mapping output         /mapping/result/
annotations            /annotations/
annotation input       /annotations/input/
annotation output      /annotation/result/
error logs             /franklin_errors/
====================== ========================

Sequence files should be fasta or sanger fastq, fasta for the ones without quality and fastq for the ones with quality. The file extension should reflect the sequence file format, fasta and sfastq.

The sequence file names for the reads should define several tags: lb (library), sm (sample) and pl (platform/technology). This convention is required if a mapping analysis is to be done with this sequence files. These tags follow the sam file header `specification <http://samtools.sourceforge.net/SAM1.pdf>`_. An example of some valid sequence file names is::

  $ ls reads/original/
  lb_mos.pl_illumina.sm_mos.sfastq  lb_pep.pl_illumina.sm_pep.sfastq
  lb_mu16.pl_454.sm_mu16.sfastq     lb_upv196.pl_454.sm_upv196.sfastq


.. _clean-reads:

Cleaning sequence reads
-----------------------

Introduction
____________

franklin can clean sanger, 454 and illumina sequences. This process usually involves vector and adaptor removal, bad quality regions trimming and short sequence filtering. There are three cleaning pipelines defined in franklin that are used depending on the platform and on the quality availability:

long reads with quality
  for sanger and 454 sequences with quality information

long reads without quality
  for sanger reads without quality information

solexa
  for short illumina reads

A collection of cleaning steps are available that compose each one of these pipelines. These steps are:

adaptor removal
  Each sequence is align against the adaptors found in a fasta file. The external tool used to do the matching is exonerate. If a match is found this section of the read is removed.

precise vector removal
  If the vector and cloning site is known lucy can be used to remove the vector in a precise way.

bad quality trimming
  There are two algorithms used to remove the bad quality sequence extremes. If the sequence is long lucy (454 and sanger) is used for this task otherwise franklin does the job (illumina).

general vector removal
  The reads are compared against the Univec database using blast to look for remaining vectors.

low complexity masking
  The regions with a low complexity are masked by using mdust

word removal
  If we know that some particular sequence appear at the beginning of the read we can remove it using this module

edge removal
  After all the other modules are run we can delete a fixed amount of bases from the sequence extremes

short sequence filtering
  When the process for one sequence is completed a minimum length criteria is applied.

The pipelines are:

long reads with quality
  adaptor removal, precise vector removal, bad quality trimming, general vector removal, low complexity masking, word removal, edge removal, and short sequence filtering 

long reads without quality
  general vector removal, bad quality trimming, low complexity masking, word removal, edge removal, and short sequence filtering 

solexa
  adaptor removal, bad quality trimming,  and short sequence filtering

Input and output files
______________________

The reads to be cleaned should be in the project directory under /reads/original/. The `naming conventions`_ should be followed by these files, especially the bit regarding to the extension. The output files will have the same names, but they will be located at /reads/cleaned/. The analysis will proceed for all sequence files found in /reads/original, if a matching file is not found in /reads/cleaned/ a new cleaned file will be generated. If a matching file is found in /reads/cleaned/ these file will not be overwritten, so the analysis for this file will not be repeated until the file from /reads/cleaned is removed.


Configuration parameters
________________________

The configuration for the cleaning analysis is found in the Cleaning section on the franklin.conf file. The parameters are:

vector_database
  The blast database that will be used to look for clonning vectors.

adaptors_file_454
  A path to a fasta file containing the adaptors used to build the 454 library. They will be removed from the cleanend reads.

adaptors_file_sanger
  Idem for the sanger sequences

adaptors_file_illumina
  Idem for the illumina sequences

words_to_remove_454
  A list of words to be removed if they are found at the start of the 454 sequences.
  
words_to_remove_sanger
  Idem for the sanger sequences

words_to_remove_illumina
  Idem for the illumina sequences
 
edge_removal -> 454_left
  A fixed number of bases to be removed from the left edge of the 454 reads.

edge_removal -> 454_right
  Idem for the right edge of the 454 reads

edge_removal -> sanger_left
  Idem for the left edge of the sanger reads

edge_removal -> sanger_right
  Idem for the right edge of the sanger reads

edge_removal -> illumina_left
  Idem for the left edge of the illumina reads

edge_removal -> illumina_right
  Idem for the right edge of the illumina reads

lucy_settings
  A path to a lucy settings file with the splice and vectors files to be used by lucy

lucy settings
_____________

The lucy settings file should have the following format::

  {'library1':{'vector_file':'lib1_vector.fasta', 'splice_file':'lib1_splice.fasta'},
   'library2':{'vector_file':'lib2_vector.fasta', 'splice_file':'lib2_splice.fasta'},}

In this file the paths to the vector and splice files for lucy should be stated for every library to be cleaned by lucy. The library name will be scraped from the read sequence file (that should follow the `naming conventions`_. The vector file is just a fasta file, the information to be set in the splice file should is explained in the lucy man page.


.. _mira-assembly:

Mira assembly
-------------

Introduction
____________

The `mira <http://sourceforge.net/apps/mediawiki/mira-assembler/index.php?title=Main_Page>`_ assembler is used to create a set of contigs with the sequence reads. The reads can come from 454, sanger and illumina sequencing. Hybrid assembly are possible. For mira configuration details refer to its documentation.

Input and output files
______________________

The input files required to do a mira analysis are the reads located in reads/cleaned. The reads files should follow the `naming conventions`_. 

Configuration parameters
________________________

The default configuration is tailored to EST assemblies. To modify the mira command line parameters you should go to the mira section in the franklin.conf file. The options are:

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

The mira assembly analysis is divided into three franklin analyses: prepare_mira_assembly, mira_assembly and select_last_assembly.

The analysis prepare_mira_assembly will create the files required as input by mira in the directory assembly/input/. These files will be created taking the reads from reads/cleaned/.

The mira_assembly analysis runs mira and creates the contigs. The files created by this analysis will be located at a timestamped directory located in assembly/. Several assemblies could be created with different parameters and each one would go into a different timestamped directory. Inside these directories a result subdirectory is created with the relevant result files.

The select_last_assembly will just make a soft link named assembly/result that points to the result subdirectory located in the latest timestamped assembly.

