
File formats
============

There are a lot of file sequence file formats. The include different information about the sequence and have a very different aspect. The most common file formats in the NGS world are: sff, fastq and fasta.

sff
---

The SFF (Standard Flowgram Format) files are the 454 equivalent to the ABI chromatogram files. The hold the information about the flowgram, the called sequence and the recommended quality clipping. These are binary files. We can obtain a more usable sequence text file by using the tools provided by Roche with the 454 machine. Alternatively we can use the `sff_extract <http://bioinf.comav.upv.es/sff_extract/index.html>` tool to obtain a fasta file.


sanger fastq
------------

The fastq format was developed to provided a convenient way of storing the sequence and the quality in the same file. These are text files and the look like::

  @seq_1
  GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
  +
  !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
  @seq_2
  ATCGTAGTCTAGTCTATGCTAGTGCGATGCTAGTGCTAGTCGTATGCATGGCTATGTGTG
  +
  208DA8308AD8SF83FH0SD8F08APFIDJFN34JW830UDS8UFDSADPFIJ3N8DAA

In this file every sequence has 4 lines. In the first line we get the name and, optionally, the description after the symbol "@". The second line has the sequence and the fourth line has the quality scores encoded as letters.


illummina fastq
---------------

This file is almost identical to a sanger fastq file, but the encoding for the quality scores is different. When we deal with a fastq file we have to be sure about which kind of file we are dealing with, an illumina fastq or a sanger fastq. Unfortunately they are not easy to distinguish. Also you have to take into account that solexa used to had a third fastq format, the solexa fastq.


