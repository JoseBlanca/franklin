.. franklin documentation master file, created by
   sphinx-quickstart on Fri Apr  9 09:16:18 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


franklin
======================

franklin is a bioinformatic application created to work on sequence analysis by using NGS (Next Generation Sequencing)  and sanger sequences. It is capable of cleaning reads, do de novo assembly or mapping against a reference and annotate SNPs, SSRs, ORFs, GO terms and sequence descriptions.
Our laboratory is focused on transcriptomic analysis, so the tool has been used and tested on transcriptomes. Some analyses will be useful for genome analysis, but since our work deals mainly with transcriptomes design tradeoffs in franklin reflect this. 

For the analyses, in most cases, franklin uses external software like: mira_, bwa_, samtools_, picard_, etc.

The application works on Linux so experience with that operating system is advisable when using it.

.. toctree::
   :maxdepth: 2
   :hidden:

   introduction
   using_franklin
   annotation
   architecture


.. _mira: http://sourceforge.net/apps/mediawiki/mira-assembler
.. _bwa: http://bio-bwa.sourceforge.net/
.. _samtools: http://samtools.sourceforge.net/
.. _picard: http://picard.sourceforge.net/index.shtml
