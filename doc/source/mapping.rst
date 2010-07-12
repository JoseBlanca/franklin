
Mapping
-------

A set of read files can be mapped against a reference genome. For the mapping ngs_backbone uses bwa with two algorithms, one for the long reads (sanger and 454) and other for the short reads (illumina). The result is a set of bam files one for each input read file or a merged bam file with all reads in it.

Input and output files
______________________

The read files should be located in reads/cleaned/ and should follow the :doc:`naming conventions <introduction>`. It is very important to set in the read file names the library, sample and platforms, otherwise the realignment and the SNP calling will fail. The reference genome should be located in mapping/reference as a fasta file.

Once bwa is finished a timestamped mapping directory will contain a bams/by_readgroup subdirectory with one bam file for each input read file. Every one of such bam files is considered to be a read group. After the mapping is finished a merge_bam analysis can be done. That analysis will merged all bam files located in bams/by_reagroup and will create an unique bam file in bams/merged.bam. This bam file will contain as many read groups as bam files are merged. Every read group will retain the information about the library, sample and platform.

Running the analysis
____________________

The analysis is run in divided in two ngs_backbone analysis:

mapping
  It maps the reads with bwa creating one bam for every input file

merge_bam
  It merges all bam files located in mapping/bams/by_readgroup into mapping/bams/merged.bam. The obtained bam will comply not only with the samtools standard but also with the picard and GATK requirements.


.. _bam-realignment:

Bam realignment
---------------

This analysis does a GATK realignment_. The mappings are usually done aligning each read with the reference genome at a time. These methodology can cause artifacts in the multiple sequence alignment obtained. GATK is capable of solving these artifacts. Their algorithm is described in its own site.

Input and output files
______________________

The only one input file should be mapping/bams/merged.bam. This bam file contains all the reads mapped to the reference genome. The output file will be also mapping/bams/merged.bam (a versioned copy).

Running the analysis
____________________

The corresponding ngs_backbone is realign_bam.


.. include:: links.txt

