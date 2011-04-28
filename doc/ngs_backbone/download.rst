
Download
========

The recommended method to get ngs_backbone is to download one of the releases, but if you want to access to the latest code you can do it at github_.

Virtual Machine
---------------

To ease the evaluation of ngs_backbone without having to follow the complete install process we have prepared a `VirtualBox <http://www.virtualbox.org/>`_ `virtual machine <http://bioinf.comav.upv.es/_downloads/ngs_machine_v2.tar.gz>`_ already pre-installed with the software and with the files required to do our `NGS workshop <http://bioinf.comav.upv.es/courses/ngs_workshop/>`_.
Be aware that this machine won't be able to run any real sized experiment and it will only be useful for testing and evaluation purposes.

ngs_backbone 1.3.1
------------------

:download:`ngs_backbone 1.3.1 <downloads/ngs_backbone-1.3.1.tar.gz>`. Relased on 28-04-2011.

Bug fix release.

ngs_backbone 1.3.0
------------------

:download:`ngs_backbone 1.3.0 <downloads/ngs_backbone-1.3.0.tar.gz>`. Relased on 15-04-2011.

Changes:

  * The vector cleaning algorithms have been reworked, now we use blast short sequence algorithm.
  * Added quality trimming for Solid reads.
  * MAF threshold added to SNP calling.
  * Improved SAM statistics.
  * Third party tools are now included in ngs_backbone to ease installation.
  * Several internal refactorings and bug fixes.


ngs_backbone 1.2.0
------------------

:download:`ngs_backbone 1.2.0 <downloads/ngs_backbone-1.2.0.tar.gz>`. Relased on 04-11-2010.

Changes:

 * Solid support added
 * Seqio csfasta support
 * Quality filtering
 * Solid mapping with bwa
 * Improved microsatellite statistics
 * Several bug fixes.
 * Experimental support for SNP effect on proteins.


ngs_backbone 1.1.0
------------------

:download:`ngs_backbone 1.1.0 <downloads/ngs_backbone-1.1.0.tar.gz>`. Relased on 31-08-2010.

Changes:

 * More stringent SNV caller for low quality sequences.
 * Clearer error messages for some common misconfigurations.
 * Several bug fixes.

ngs_backbone 1.0.0
------------------

:download:`ngs_backbone 1.0.0 <downloads/ngs_backbone-1.0.0.tar.gz>`. Relased on 12-07-2010.

Changes:

 * box plots added for for the sequence quality along the sequence.
 * Several bug fixes.

ngs_backbone 0.3
----------------

:download:`ngs_backbone 0.3.0 <downloads/ngs_backbone-0.3.0.tar.gz>`. Relased on 21-06-2010.

Changes:

 * Annotation statistics analysis added.
 * The CAP enzymes are now written to the vcf file.
 * SNV structure simplified to save memory (incompatible change with old projects).
 * Database sequences moved from repr to pickle (memory improvement).
 * Improved setup with dependency check.
 * Install documentation expanded.
 * Several bug fixes.

ngs_backbone 0.2
-----------------

:download:`ngs_backbone 0.2.0 <downloads/ngs_backbone-0.2.0.tar.gz>`. Relased on 3-06-2010.

Changes:

 * ngs_backbone is now a parallel application.

 * Now we look for the short adaptors by using regular expressions.

 * bam to sam conversion is now much faster.

 * Statistics for the bam files added.

 * documentation for an NGS workshop using ngs_backbone added.

 * UniVec use is now optional.

 * Some directory names were changed.

 * Better interoperability with IGV.

 * New SNV filter added to filter out SNVs not read in enough samples or libraries.

 * Restriction enzyme search moved from EMBOSS to Biopython to improve speed.

 * Statistics module refactored to improve speed.

 * GO annotation can can work with the annotation file skipping the b2gpipe.

 * Config module refactored to work with old ConfigObj versions.

 * numpy dependency has been removed.

 * The soft links are now relative.

 * The select last analysis are now run automatically.

 * Other minor performance improvements.

 * Lots of bug fixes.

ngs_backbone 0.1
-----------------

:download:`ngs_backbone 0.1.0 <downloads/ngs_backbone-0.1.0.tar.gz>`. Relased on 30-04-2010.

Changes:

 * Initial Release.

Feedback, positive or negative will be appreciated.


.. include:: links.txt

