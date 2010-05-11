
seq_io
======

seq_io.py is a little utility distributed with ngs_backbone that allows us to move sequence files from one sequence format to other. To change between formats we just have to tell seq_io.py which are the input and output files and which are the input and output formats. The parameters to use are::

  $ seqio.py -h
  Usage: seqio.py [options]

  Options:
    -h, --help            show this help message and exit
    -s INSEQFILE, --inseqfile=INSEQFILE input sequence file
    -q INQUALFILE, --inqualfile=INQUALFILE input quality file
    -f INFORMAT, --informat=INFORMAT input file format
    -o OUTSEQFILE, --outseqfile=OUTSEQFILE output sequence file
    -l OUTQUALFILE, --outqualfile=OUTQUALFILE output quality file
    -t OUTFORMAT, --outformat=OUTFORMAT output file format

For instance, to go from sanger fastq to fasta and qual we would do::

  seqio.py -s seq.sfastq -f sfastq -o seq.fasta -l seq.qual -t fasta

And to do the reverse::

  $ seqio.py -s seq.fasta -q seq.qual -f fasta -o seq.sfastq -t sfastq

