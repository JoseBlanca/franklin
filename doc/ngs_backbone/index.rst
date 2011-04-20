.. ngs_backbone documentation master file, created by
   sphinx-quickstart on Fri Apr  9 09:16:18 2010.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


ngs_backbone
============

ngs_backbone is a bioinformatic pipeline created to work on Next Generation Sequence (NGS) analysis as well as with sanger sequences.
It is capable of cleaning reads, prepare a novo assembly, map reads against a reference, look for SNPs and SSRs, and do some function annotation like :ORFs, GO terms and sequence descriptions.

ngs_backbone:

  * is mainly a wrapper around external tools like: bwa_, samtools_, blast, etc.
  * uses standard files like fastq_, SAM_, VCF_ and GFF_ to ease the interoperability with other tools.
  * is a command line application.
  * is free software released under the AGPL license with the hope that could be useful to other laboratories.
  * runs in Linux so experience with that operating system is advised when using it.
  * is written in Python.
  * can run in parallel using single node multicore systems and computer clusters.

If you find any problem when running ngs_backbone we'd love to hear from you, so please inform us.


Help
----

To get help there is a mailing list for ngs_backbone.
If you have any question send an email to the list: ngs_backbone@upv.es

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

.. _fastq: http://en.wikipedia.org/wiki/FASTQ_format
.. _bwa: http://bio-bwa.sourceforge.net/
.. _samtools: http://samtools.sourceforge.net/
.. _SAM: http://samtools.sourceforge.net/SAM-1.3.pdf
.. _VCF: http://www.1000genomes.org/wiki/Analysis/vcf4.0
.. _GFF: http://www.sequenceontology.org/resources/gff3.html

