
Installation
============

This is a summary on how to install ngs_backbone, below you have a :ref:`detailed explanation <step_by_step>`.

If you just want to check out how ngs_backbone works you could also download a complete `VirtualBox <http://www.virtualbox.org/>`_ `virtual machine <http://bioinf.comav.upv.es/_downloads/ngs_machine_v3.tar.gz>`_ already pre-installed with the software and with the files required to do our `NGS workshop <http://bioinf.comav.upv.es/courses/sequence_analysis/index.html>`_.
Although this machine will not be powerful enough to run most real sized analyses it can be an easy way to get a felling of the analysis process without spending any time at all setting up the environment.

To install ngs_backbond python 2.6 or better is required.
Also you need the python libraries Biopython_, and ConfigObj_.
An easy way to install these libraries would be to use the packaging system of your Linux distribution or alternatively just run the following commands::

  $ easy_instal biopython
  $ easy_install configobj

If you don't have easy_install installed in your computer.
You can grab it from here: setuptools_.
Other optional, but highly recommended, dependencies are: psubprocess_, pysam_ and matplotlib_.
All of these libraries are installed in the same way::

  $ tar -zxvf package.tar.gz
  $ cd package
  $ python setup.py install
 
Most of the third party tools required by ngs_backbone are included in it, but if you're planing to annotate GO terms you will have to install b2g4pipe(blast2go_).

For the vector cleaning the Univec_ database can be used.
If you configure the cleaning process to use it you also have to install it.

.. _step_by_step:

Step by step installation instructions
--------------------------------------

python tools
____________

ngs_backbone requires python2.6.
Chance are that it is already included in your distribution, but if it is not your distribution you can `download <http://www.python.org/download/releases/>`_ it and install it.
The other requirements are Python libraries.
Biopython and ConfigObj are required and pysam, psubprocess and matplotlib are optional.

It is quite possible that your distribution already includes Biopython and ConfigObj so you can use your package manager to install them, but we are going to explain here the manual process.
The install process is simple, you just have to download a bunch of python tools, unpack them and run "python2.6 install" on them.

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

To call the SNP you will need the library pysam_.

::

  $ tar -xvzf pysam-0.4.1.tar.gz
  $ cd pysam-0.4.1
  $ python2.6 setup.py install 

To create the charts for the statistics you will need matplotlib_.

::

  $ tar -xvzf matplotlib-1.0.1.tar.gz
  $ cd matplotlib-1.0.1
  $ python2.6 setup.py install

Once we have it all we can install ngs_backbone.

::

  $ tar -xvzf ngs_backbone.0.2.0.tar.gz
  $ cd ngs_backbone.0.2.0
  $ python2.6 setup.py install


Java tools
__________

Three java tools are used: picard_ ,GATK_ and blast2go_.
The first two are provided with ngs_backbone so you don't need to worry about them.

To install blast2go_, just download it and unpack it.

::

  $ tar zxvf b2g4pipe235.tar.gz
  $ mv b2g4pipe235 /usr/local
  $ updatedb

It is advisable to run updatedb after setting everything to ease the ngs_backbone configuration.

Finally to run the java tools you need java.
If you don't have it already installed just do.

::

  $ apt-get install sun-java6-jre


After installing the whole pipeline you can run the `NGS workshop tutorial <http://bioinf.comav.upv.es/courses/ngs_workshop/>`_ to test the whole system.


.. include:: links.txt

