

.. _clean-reads:

Cleaning sequence reads
-----------------------

ngs_backbone can clean sanger, 454 and illumina sequences. This process usually involves vector and adaptor removal, bad quality regions trimming and short sequence filtering. There are three cleaning pipelines defined in ngs_backbone that are used depending on the platform and on the quality availability:

long reads with quality
  for sanger and 454 sequences with quality information

long reads without quality
  for sanger reads without quality information

solexa
  for short illumina reads

A collection of cleaning steps are available that compose each one of these pipelines. These steps are:

adaptor removal
  Each sequence is align against the adaptors found in a fasta file. The external tool used to do the matching is exonerate. If a match is found this section of the read is removed. Short adaptors will be treated as such.

short adaptor removal
  ngs_backbone will look for adaptors shorter than 15 bp with exact matches.

precise vector removal
  If the vector and cloning site is known lucy can be used to remove the vector in a precise way.

bad quality trimming
  There are two algorithms used to remove the bad quality sequence extremes. If the sequence is long lucy (454 and sanger) is used for this task otherwise ngs_backbone does the job (illumina).

bad quality trimming by Ns
  When there are no qualities available we can infer the quality from the percent of Ns found in the sequence. This step removes regions with a high number of Ns.

general vector removal
  The reads are compared against the Univec database using blast to look for remaining vectors.

low complexity masking
  The regions with a low complexity are masked by using mdust

edge removal
  After all the other modules are run we can delete a fixed amount of bases from the sequence extremes

short sequence filtering
  When the process for one sequence is completed a minimum length criteria is applied.

The pipelines are:

long reads with quality
  adaptor removal, precise vector removal, bad quality trimming, general vector removal, low complexity masking, word removal, edge removal, and short sequence filtering

long reads without quality
  general vector removal, bad quality trimming by Ns, low complexity masking, word removal, edge removal, and short sequence filtering

solexa
  adaptor removal, bad quality trimming,  and short sequence filtering

Input and output files
______________________

The reads to be cleaned should be in the project directory under /reads/raw/. The :doc:`naming conventions <introduction>` should be followed by these files, especially the bit regarding to the extension. The output files will have the same names, but they will be located at /reads/cleaned/. The analysis will proceed for all sequence files found in /reads/raw, if a matching file is not found in /reads/cleaned/ a new cleaned file will be generated. If a matching file is found in /reads/cleaned/ these file will not be overwritten, so the analysis for this file will not be repeated until the file from /reads/cleaned is removed.

.. _clean-config:

Configuration parameters
________________________

The configuration for the cleaning analysis is found in the Cleaning section on the ngs_backbone.conf file. The parameters are:

vector_database
  The blast database that will be used to look for clonning vectors.

adaptors_file_454
  A path to a fasta file containing the adaptors used to build the 454 library. They will be removed from the cleanend reads.

adaptors_file_sanger
  Idem for the sanger sequences

adaptors_file_illumina
  Idem for the illumina sequences

short_adaptors_454
  A list of words to be removed. They can be regular expressions.

short_adaptors_sanger
  Idem for the sanger sequences

short_adaptors_illumina
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

strip_n_percent
  Threshold used for the trimming of regions with a lot of Ns in sequences with no quality information available. Lower values (e.g. 1.5) are more stringent.

min_seq_length
  The minimum sequence length allowable after the cleaning is done. All sequences shorter than these values will be discarded. This is a subsection with one value for each platform 454, sanger and illumina.

lucy -> vector_settings
  A path to a lucy settings file with the splice and vectors files to be used by lucy

lucy ->bracket
    Look at lucy man page before changing defaults.

lucy -> window
    Look at lucy man page before changing defaults.

lucy -> error
    Look at lucy man page before changing defaults.


lucy settings
_____________

The lucy settings file should have the following format:

::

  {'library1':{'vector_file':'lib1_vector.fasta', 'splice_file':'lib1_splice.fasta'},
   'library2':{'vector_file':'lib2_vector.fasta', 'splice_file':'lib2_splice.fasta'},}

In this file the paths to the vector and splice files for lucy should be stated for every library to be cleaned by lucy. The library name will be scraped from the read sequence file (that should follow the :doc:`naming conventions <introduction>`. The vector file is just a fasta file, the information to be set in the splice file should is explained in the lucy man page.


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

Mapping
-------

A set of read files can be mapped against a reference genome. For the mapping ngs_backbone uses bwa with two algorithms, one for the long reads (sanger and 454) and other for the short reads (illumina). The result is a set of bam files one for each input read file or a merged bam file with all reads in it.

Input and output files
______________________

The read files should be located in reads/cleaned/ and should follow the :doc:`naming conventions <introduction>`. It is very important to set in the read file names the library, sample and platforms, otherwise the realignment and the SNP calling will fail. The reference genome should be located in mapping/reference as a fasta file.

Once bwa is finished a timestamped mapping directory will contain a result/by_readgroup subdirectory with one bam file for each input read file. Every one of such bam files is considered to be a read group. After the mapping is finished a merge_bam analysis can be done. That analysis will merged all bam files located in result/by_reagroup and will create an unique bam file in result/merged_bam. This bam file will contain as many read groups as bam files are merged. Every read group will retain the information about the library, sample and platform.

Running the analysis
____________________

The analysis is run in divided in three ngs_backbone analysis:

mapping
  It maps the reads with bwa creating one bam for every input file

select_last_mapping
  It creates a soft link from mapping/result to mapping/last_timestamped_mapping/result

merge_bam
  It merges all bam files located in mapping/result/by_readgroup into mapping/result/merged.bam. The obtained bam will comply not only with the samtools standard but also with the picard and GATK requirements.


.. _bam-realignment:

Bam realignment
---------------

This analysis does a `GATK <http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit>`_ `realignment <http://www.broadinstitute.org/gsa/wiki/index.php/Local_realignment_around_indels>`_. The mappings are usually done aligning each read with the reference genome at a time. These methodology can cause artifacts in the multiple sequence alignment obtained. GATK is capable of solving these artifacts. Their algorithm is described in its own site.

Input and output files
______________________

The only one input file should be mapping/result/merged.bam. This bam file contains all the reads mapped to the reference genome. The output file will be also mapping/result/merged.bam (a versioned copy).

Running the analysis
____________________

The corresponding ngs_backbone is realign_bam.



