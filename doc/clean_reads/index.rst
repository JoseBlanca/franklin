.. clean_reads documentation master file, created by
   sphinx-quickstart on Wed Apr 20 16:03:35 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

clean_reads
===========

clean_reads cleans NGS (Sanger, 454, Illumina and solid) reads.
It can trim:

  * bad quality regions,
  * adaptors,
  * vectors, and
  * regular expresssions.

It also filters out the reads that do not meet a minimum quality criteria based on the sequence length and the mean quality.

It uses several algorithms and third party tools to carry out the cleaning.
The third party tools used are: lucy_, blast_, mdust_ and trimpoly_.

The functionality offered by clean_reads is similar to the cleaning capabilities of the ngs_backbone_ pipeline.
In fact, both tools use the same code base and are just different interfaces on top of a Python library called franklin_.

*Help*

To get help you can use the ngs_backbone maling list.
If you have any question send an email to the list: ngs_backbone@upv.es

To subscribe, unsubscribe or see the archive `list <https://listas.upv.es/mailman/listinfo/ngs_backbone>`_.

.. include:: ../ngs_backbone/links.txt

.. toctree::
   :maxdepth: 2
   :hidden:

   install
   manual
   download


