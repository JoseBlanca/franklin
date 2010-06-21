.. ngs_backbone documentation master file, created by
   sphinx-quickstart on Fri Apr  9 09:16:18 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


ngs_backbone
============

ngs_backbone is a bioinformatic application created to work on sequence analysis by using NGS (Next Generation Sequencing)  and sanger sequences. It is capable of cleaning reads, do de novo assembly or mapping against a reference and annotate SNPs, SSRs, ORFs, GO terms and sequence descriptions.
Our laboratory is focused on transcriptomic analysis, so the tool has been used and tested on transcriptomes. Some analyses will be useful for genome analysis, but since our work deals mainly with transcriptomes design tradeoffs in ngs_backbone reflect this.

ngs_backbone can run in parallel using single node multicore systems and computer clusters.

For the analyses, in most cases, ngs_backbone uses external software like: mira_, bwa_, samtools_, picard_, etc.

The application works on Linux so experience with that operating system is required when using it.

ngs_backbone is free software so contributions and shared development will be welcomed.

ngs_backbone is at this time beta software so it will not be bug free. It works for us and we hope that it will work for you, but if you find any problem, please inform us.

.. toctree::
   :maxdepth: 2
   :hidden:

   introduction
   install
   cleaning
   assembly
   mapping
   annotation
   snv_filters
   tutorials
   ngs_workshop <ngs_workshop/index>
   seqio
   architecture
   download


.. _mira: http://sourceforge.net/apps/mediawiki/mira-assembler
.. _bwa: http://bio-bwa.sourceforge.net/
.. _samtools: http://samtools.sourceforge.net/
.. _picard: http://picard.sourceforge.net/index.shtml
