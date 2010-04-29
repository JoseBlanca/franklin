
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
    -t OUTSEQFILE, --outseqfile=OUTSEQFILE output sequence file
    -r OUTQUALFILE, --outqualfile=OUTQUALFILE output quality file
    -e OUTFORMAT, --outformat=OUTFORMAT output file format

For instance, to go from sanger fastq to fasta and qual we would do::

  seqio.py -s seq.sfastq -f sfastq -t seq.fasta -r seq.qual -e fasta

And to do the reverse::

  $ seqio.py -s seq.fasta -q seq.qual -f fasta -t seq.sfastq -e sfastq

