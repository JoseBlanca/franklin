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

If you find any problem when running ngs_backbone please inform us.


**Help**

To get help there is a mailing list to help with ngs_backbone. If you have any question send an email to the list: ngs_backbone@upv.es

To subscribe, unsubscribe or see the archive `list <https://listas.upv.es/mailman/listinfo/ngs_backbone>`_.


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
   seqio
   architecture
   download


.. _mira: http://sourceforge.net/apps/mediawiki/mira-assembler
.. _bwa: http://bio-bwa.sourceforge.net/
.. _samtools: http://samtools.sourceforge.net/
.. _picard: http://picard.sourceforge.net/index.shtml
