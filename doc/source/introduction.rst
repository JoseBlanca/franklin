
Introduction
============

ngs_backbone is an easy to use application capable of doing in a reliable way some NGS analyses. The main design directives have been:

 * analyses reproducibility
 * ease of use
 * modularity
 * standard format output

By using ngs_backbone we can run a complete analysis in a reproducible way. Every analysis parameter is configured in a text file and a log file is generated as the analysis progresses. Also the analysis is easy to do, instead of running multiple scripts and programs to do an analysis only one command is required.

An application like that in a fast moving field, as the NGS is, has the risk of stagnate rapidly. To avoid this pitfall special care has been taken in the design of the architecture of ngs_backbone. Everything within the application is modular and several independent layers have been differentiated within the code to facilitate the maintenance of the code. In fact before its public release ngs_backbone has gone through refactorings of several modules that have not affected the overall structure of the application.

Usage
=====

The main ideas to consider when using ngs_backbone are: project and analysis. A *project* is a directory (with its subdirectories) that includes a configuration file and all input and output files. An *analysis* takes some inputs from the project and creates some outputs. ngs_backbone knows where the input and output files are for every analysis because the project directory structure is the same for every project. For instance the reads to clean are always in the directory /reads/original and the cleaned reads always are in /reads/cleaned/.

The configuration parameters required for every analysis are stored in the configuration file (ngs_backbone.conf) located in the project directory. A sample of this configuration file is created by ngs_backbone when a new project is created. This file should be tweaked to adapt the analyses to your requirements before running them.

ngs_backbone has only two main executables: backbone_create_project.py and backbone_analysis.py. The first one is used to create a project from scratch and the second one to run all analyses.

Creating a new project
----------------------

::

  $ backbone_create_project.py -p project_name
  $ ls -l project_name
  total 4
  -rw-r--r-- 1 jose jose 2479 abr 16 09:30 ngs_backbone.conf

This command will create a new directory named project_name with a file named ngs_backbone.conf in it. These are the two hallmarks that define a ngs_backbone *project*, the directory and the configuration file.

Running an analysis
-------------------

::

  $ backbone_analysis.py -a analysis_name

This command will run the analysis on the data present in the project using the parameters found in the configuration file. The output files will also be located in the project.

.. _naming:

Naming conventions
==================

The ngs_backbone usage is heavily based on directory and file name conventions. If the input files are not located where ngs_backbone expects to find them or they have non-standard file names the analysis will fail.

====================== ========================
file or directory type location
====================== ========================
configuration file     ngs_backbone.conf
log file               ngs_backbone.log
raw reads              reads/raw/
clean reads            reads/cleaned/
assemblies             assembly/
assembly input         assembly/input/
assembly output        assembly/result/
mappings               mapping/
mapping output         mapping/bams/
annotations            annotations/
annotation input       annotations/input/
annotation output      annotation/features/
error logs             backbone_errors/
====================== ========================

Sequence files should be fasta or sanger fastq, fasta for the ones without quality and fastq for the ones with quality. The file extension should reflect the sequence file format, fasta and sfastq.

The sequence file names for the reads should define several tags: lb (library), sm (sample) and pl (platform/technology). This convention is required if a mapping analysis is to be done with this sequence files. These tags follow the sam file header `specification <http://samtools.sourceforge.net/SAM1.pdf>`_. An example of some valid sequence file names is::

  $ ls reads/raw/
  lb_mos.pl_illumina.sm_mos.sfastq  lb_pep.pl_illumina.sm_pep.sfastq
  lb_mu16.pl_454.sm_mu16.sfastq     lb_upv196.pl_454.sm_upv196.sfastq

Available analyses
==================

The available analyses are:

========================================    =================================================
analysis                                    description
========================================    =================================================
:ref:`clean-reads`                          sequence reads cleaning
:ref:`mira-assembly`                        Assembly reads into a contig set with  mira
:ref:`mapping`                              bwa read mapping against a reference sequence
:ref:`bam-realignment`                      GATK bam realignment
:ref:`snp-calling`                          SNP annotation from a bam file
:ref:`orf-annotation`                       ESTScan ORF annotation
:ref:`ortholog-annotation`                  reciprocal blast based ortholog annotation
:ref:`description-annotation`               description blast based annotation
:ref:`ssr-annotation`                       microsatellite sputnik based annotation
:ref:`intron-annotation`                    est2genome cDNA based annotation
:ref:`go-annotation`                        blast2go  annotation
========================================    =================================================

Parallel operation
==================

Running ngs_backbone with multiple subprocesses is as easy as setting the configuration option threads to True. ngs_backbone will run will as many subprocesses as cpu cores are found in the computer. Also the threads option can be set to an integer and ngs_backbone will run with as many subprocess as indicated.



