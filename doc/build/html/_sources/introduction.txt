
Description
===========

franklin is a bioinformatic application created to work on sequence analysis by using NGS (Next Generation Sequencing)  and sanger sequences. It is capable of cleaning reads, do de novo assembly or mapping against a reference and annotate SNPs, SSRs, ORFs, GO terms and sequence descriptions.
Our laboratory is focused on transcriptomic analysis, so the tool has been used and tested on transcriptomes. Some analyses will be useful for genome analysis, but since our work deals mainly with transcriptomes design tradeoffs in franklin reflect this. 

For the analyses, in most cases, franklin uses external software like: mira_, bwa_, samtools_, picard_, etc.

The application works on Linux so experience with that operating system is advisable when using it.

Philosophy
==========

franklin is an easy to use application capable of doing in a reliable way some NGS analyses. The main design directives have been:

 * analyses reproducibility
 * ease of use
 * modularity
 * standard format output

By using franklin we can run a complete analysis in a reproducible way. Every analysis parameter is configured in a text file and a log file is generated as the analysis progresses. Also the analysis is easy to do, instead of running multiple scripts and programs to do an analysis only one command is required.

An application like that in a fast moving field, as the NGS is, has the risk of stagnate rapidly. To avoid this pitfall special care has been taken in the design of the architecture of franklin. Everything within the application is modular and several independent layers have been differentiated within the code to facilitate the maintenance of the code. In fact before its public release franklin has gone through refactorings of several modules that have not affected the overall structure of the application.



Installation
============

To install franklin python 2.6 and Biopython_ are required. python is installed by default in the usual Linux distributions but you should check the version. Installing Biopython_ is as easy as installing franklin. Once the tarball is downloaded you have to run the following command::

  $ python setup.py install

franklin requires several external tools to run the analyses. So before running the analyses please install them.

=============================  ================================================================
analysis                       external tools required
=============================  ================================================================
:ref:`clean-reads`             lucy_, exonerate_, blast_, Univec_ (database), mdust_, trimpoly_
:ref:`mira-assembly`           mira_
:ref:`mapping`                 bwa_, samtools_, picard_
:ref:`bam-realignment`         GATK_
:ref:`snp-calling`             pysam_
:ref:`orf-annotation`          ESTScan_
:ref:`ortholog-annotation`     blast_
:ref:`description-annotation`  blast_
:ref:`ssr-annotation`          sputnik
:ref:`intron-annotation`       blast_, emboss_
:ref:`go-annotation`           blast_, blast2go_
=============================  ================================================================



.. _mira: http://sourceforge.net/apps/mediawiki/mira-assembler
.. _bwa: http://bio-bwa.sourceforge.net/
.. _samtools: http://samtools.sourceforge.net/
.. _picard: http://picard.sourceforge.net/index.shtml
.. _pysam: http://code.google.com/p/pysam/
.. _GATK: http://www.broadinstitute.org/gsa/wiki/index.php/The_Genome_Analysis_Toolkit
.. _Biopython: http://biopython.org/wiki/Main_Page
.. _lucy: http://lucy.sourceforge.net/
.. _exonerate: http://www.ebi.ac.uk/~guy/exonerate/
.. _blast: http://web.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _Univec: http://www.ncbi.nlm.nih.gov/VecScreen/UniVec.html
.. _mdust: http://compbio.dfci.harvard.edu/tgi/software/
.. _trimpoly: http://compbio.dfci.harvard.edu/tgi/software/
.. _ESTScan: http://estscan.sourceforge.net/
.. _emboss: http://emboss.sourceforge.net/
.. _blast2go: http://www.blast2go.org/

