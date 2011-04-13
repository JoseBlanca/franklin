
Installation
============

This is a summary on how to install ngs_backbone, below you have a :ref:`detailed explanation <step_by_step>`.

If you just want to check out how ngs_backbone works you could also download a complete `VirtualBox <http://www.virtualbox.org/>`_ :download:`virtual machine <downloads/ngs_machine_v2.tar.gz>` already pre-installed with the software and with the files required to do our `NGS workshop <http://bioinf.comav.upv.es/courses/ngs_workshop/>`_.
Although this machine will not be powerful enough to run most real sized analyses it can be an easy way to get a felling of the analysis process without spending any time at all setting up the environment.

To install ngs_backbond python 2.6 or better is required. Also you need the python libraries Biopython_, and ConfigObj_. You can install these libraries with the following commands:

  $ easy_instal biopython
  $ easy_install configobj

If you don't have easy_install installed in your computer, You can grab it from here: setuptools_. Other optional, but highly recommended, dependencies are: psubprocess_, pysam_ and matplotlib_. All of these libraries are installed in the same way:

  $ tar zxvf package.tar.gz
  $ cd package
  $ python setup.py install
 
python is installed by default in the usual Linux distributions but you should check the version. Installing a python library is as easy as installing ngs_backbone. Once you have the tarball downloaded just run the following command::

  $ python2.6 setup.py install

ngs_backbone requires also several external tools to run the analyses. Most of the applications are distributed in ngs_backbone so you don't need to install them. Due to license reasons we can not distribute b2g4pipe(blastogo_).

=============================  ================================================================
analysis                       external tools required
=============================  ================================================================
:ref:`clean-reads`             lucy_, exonerate_, blast_, Univec_ (database), mdust_, trimpoly_
:ref:`mira-assembly`           mira_
:doc:`mapping`                 bwa_, samtools_, picard_
:ref:`bam-realignment`         GATK_
:ref:`snp-calling`             pysam_
:ref:`orf-annotation`          ESTScan_
:ref:`ortholog-annotation`     blast_
:ref:`description-annotation`  blast_
:ref:`ssr-annotation`          sputnik
:ref:`intron-annotation`       blast_, emboss_
:ref:`go-annotation`           blast_, blast2go_
=============================  ================================================================


.. _step_by_step:

Step by step installation instructions
--------------------------------------

python tools
____________

ngs_backbone requires python2.6. If you don't have it already installed in your distribution `download <http://www.python.org/download/releases/>`_ the source code and install it. The other requirements are python libraries. Biopython and ConfigObj are required and pysam, psubprocess and matplotlib are optional.

If your distribution include python2.6 chances are that Biopython and ConfigObj might be packages by your distribution, but we are going to explain here the manual process. The install process is simple, you just have to download a bunch of python tools, unpack them and run "python2.6 install" on them.

To install Biopython you need Numpy_. Download Biopython_ and install it.

::
  $ tar -xvzf numpy-1.4.1.tar.gz
  $ cd numpy-1.4.1 

  $ tar -xvzf biopython-1.54.tar.gz
  $ cd biopython-1.54
  $ python2.6 setup.py install

Download ConfigObj_ and install it.

::

  $ unzip configobj-4.7.1.zip
  $ cd configobj-4.7.1
  $ python2.6 setup.py install

Download and install psubprocess_ (If you don't install it you won't be able to run ngs_backbone in parallel).

::

  $ tar -xvzf psubprocess.0.1.1.tar.gz
  $ cd psubprocess.0.1.1
  $ python2.6 setup.py install

To call the SNP you will need the library pysam_. pysam requires pyrex, so install it.

::

  $ tar -xvzf Pyrex-0.9.9.tar.gz
  $ cd Pyrex-0.9.9
  $ python2.6 setup.py install
  $ tar -xvzf pysam-0.2.tar.gz
  $ cd pysam-0.2
  $ python2.6 setup.py install 

To create the charts for the statistics you will need matplotlib_.

::

  $ tar -xvzf matplotlib-0.99.3.tar.gz
  $ cd matplotlib-0.99.3
  $ python2.6 setup.py install

Once we have it all we can install ngs_backbone.

::

  $ tar -xvzf ngs_backbone.0.2.0.tar.gz
  $ cd ngs_backbone.0.2.0
  $ python2.6 setup.py install


Java tools
__________

Three java tools are used to manage the sam files: picard_ ,GATK_ and blast2go_ . The two first are provided in the ngs_backbone tarball.

::

  $ apt-get install sun-java6-jre

blast2go_, just download it and unpack it.

::

  $ tar zxvf b2g4pipe235.tar.gz
  $ mv b2g4pipe235 /usr/local
  $ updatedb

It is advisable to run updatedb after setting everything to ease the ngs_backbone configuration.

After installing the whole pipeline you can run the `NGS workshop tutorial <http://bioinf.comav.upv.es/courses/ngs_workshop/>`_ to test the whole system.


.. include:: links.txt

