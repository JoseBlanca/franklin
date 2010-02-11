#!/usr/bin/env python
'''This script runs cleaning pipelines for sequence reads.

There are several pipelines to chose depending on the type of read (sanger, 454,
solexa) and on the availability of the quality data.
It takes a sequence fasta file (and optionally a quality fasta file) and it
returns the cleaned reads in a new fasta file.
The cleaning done depends on the pipeline chosen. Within the different pipelines
there are modules for:
- quality trimming based on the quality information (based on lucy or in a
  custom function),
- quality trimming based on the presence of intederminations (with trimpoly),
- low quality masking (with trimpoly),
- poly-A masking (with trimpoly)
- vector trimming using blast and a vector database like Univec,
- adaptor trimming using exonerate and a adaptor multifasta file,
- repeat masking using RepeatMasker.

The pipelines use a working directory in which several temporal files are
written. The user can specify a which working directory should be used.
If the user request it an intermediate file for every pipeline step can be
created in the given working directory. To do that the checkpoint option should
be requested.
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

import logging, os
from optparse import OptionParser
from franklin.pipelines.pipelines import seq_pipeline_runner

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -s fastafile [-q quality_file]')
    parser.add_option('-s', '--fastafile', dest='fastafile',
                      help='Sanger input fasta')
    parser.add_option('-q', '--qualfile', dest='qualfile',
                      help='Quality fasta file')
    parser.add_option('-d', '--directory', dest='directory',
                      help='Working diretory')
    parser.add_option('-v', '--vecdb', dest='vecdb',
                      help='vector database')
    parser.add_option('-a', '--adaptors', dest='adaptors',
                      help='adaptors fasta file')
    parser.add_option('-o', '--seqoutput', dest='out_seq',
                      help='sequence output file')
    parser.add_option('-u', '--qualoutput', dest='out_qual',
                      help='quality output file')
    parser.add_option('-l', '--logfile', dest='logfile',
                      help='Log file')
    parser.add_option('-p', '--pipeline', dest='pipeline',
                      help='Pipeline type. Look at the source' )
    parser.add_option('-c', '--checkpoint', action="store_true",
                      dest='checkpoint', help='Do we need checkpoints?')
    parser.add_option('-w', '--words', dest='words',
                      help='file containing words in the first column')
    parser.add_option('--lucy_vector', dest='lucy_vector',
                      help='file used by lucy for vector clipping(vector)')
    parser.add_option('--lucy_splice_site', dest='lucy_splice_site',
                     help='file used by lucy for vector clipping(splice_sites)')

    parser.add_option('-f', '--inputformat', dest="file_format",
     help='Input file format: fasta(default), fastq, fastaq, fastq-solexa, ...')

    return parser.parse_args()

def set_parameters():
    '''It sets the parameters for the script.

    For some options there are some default values that will be applied here.
    '''
    options = parse_options()[0]

    # kind
    if options.pipeline is None:
        raise RuntimeError('Pipeline is mandatory, use -p pipeline')
    else:
        pipeline = options.pipeline
    # do heckpoints?
    if options.checkpoint is None:
        checkpoint = None
    else:
        checkpoint = options.checkpoint

    ########### file format ###############################
    if options.file_format is None:
        file_format = None
    else:
        file_format = options.file_format
    ############### input output parameters ################
    io_fhands     = {}

    if options.fastafile is None:
        raise RuntimeError('Script at least needs an input file (fasta)')
    else:
        io_fhands['in_seq'] = open(options.fastafile, 'r')

    if options.qualfile is None:
        io_fhands['in_qual'] = None
    else:
        io_fhands['in_qual'] = open(options.qualfile, 'r')

    if options.out_seq is None:
        io_fhands['out_seq'] = open('result.seq.fasta', 'w')
    else:
        io_fhands['out_seq'] = open(options.out_seq, 'w')

    if options.out_qual is None:
        if  io_fhands['in_qual'] is None:
            io_fhands['out_qual'] = None
        else:
            io_fhands['out_qual'] = open('result.qual.fasta', 'w')
    else:
        io_fhands['out_qual'] = open(options.out_qual, 'w')
    ##########################################################
    # do checkpoints?
    if options.checkpoint is None:
        checkpoint = False
        work_dir   = None
    else:
        checkpoint = options.checkpoint
        if options.directory is None:
            if not os.path.exists('cleaning_pipeline'):
                os.mkdir('cleaning_pipeline')
            work_dir = 'cleaning_pipeline'
        else:
            work_dir = options.directory
            if not os.path.exists(work_dir):
                os.mkdir(work_dir)
    # logs?
    if options.logfile is None:
        log_fhand = open('clean_pipeline.log', 'w')
    else:
        log_fhand = open(options.logfile, 'w')

    ########### step configuration  ###########################
    configuration = {}
    if options.vecdb is not None:
        configuration['remove_vectors'] = {}
        configuration['remove_vectors']['vectors'] = options.vecdb
    if options.adaptors is not None:
        configuration['remove_adaptors'] = {}
        configuration['remove_adaptors']['vectors'] = options.adaptors
    if options.words is not None:
        configuration['word_masker'] = {}
        configuration['word_masker']['words'] = \
                       [line.strip().split()[0]for line in open(options.words)]
    if options.lucy_vector is not None and options.lucy_splice_site is not None:
        configuration['strip_lucy'] = {}
        configuration['strip_lucy']['vector'] = [options.lucy_vector,
                                                     options.lucy_splice_site]
    return (io_fhands, work_dir, log_fhand, pipeline, configuration, checkpoint,
            file_format)


def main():
    'The main function'

    # Set parameters
    (io_fhands, work_dir, log_fhand, pipeline, config, checkpoint,
     file_format) = set_parameters()

    #Loggin facilities
    logging.basicConfig(filename=log_fhand.name, level=logging.INFO,
                        format='%(asctime)s %(message)s')

    # Run the analisis step by step
    seq_pipeline_runner(pipeline=pipeline, configuration=config,
                        io_fhands=io_fhands, file_format=file_format)


if __name__ == '__main__':
    main()
