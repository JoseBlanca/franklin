'''
This script clean(mask and strip) and filter sequence depending in a lot of
variables in sanger sequences.
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import logging, os
from optparse import OptionParser
from biolib.seq_cleaner import pipeline_runner

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
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
                      dest='checkpoint', help='Do we need checkpoints?' )

    return parser.parse_args()

def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''
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
            if not os.path.exists('clean_sanger_tmp'):
                os.mkdir('clean_sanger_tmp')
            work_dir = 'clean_sanger_tmp'
        else:
            work_dir = options.directory
            if not os.path.exists(work_dir):
                os.mkdir(work_dir)
    # logs?
    if options.logfile is None:
        log_fhand = open('sanger_clean.log', 'w')
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

    return io_fhands, work_dir, log_fhand, pipeline, configuration, checkpoint


def main():
    'The main function'

    # Set parameters
    (io_fhands, work_dir, log_fhand, pipeline, config, checkpoint) =\
                                                        set_parameters()

    #Loggin facilities
    logging.basicConfig(filename=log_fhand.name, level=logging.INFO,
                        format='%(asctime)s %(message)s')

    # Run the analisis step by step
    pipeline_runner(pipeline, config, io_fhands, work_dir, checkpoint)


if __name__ == '__main__':
    main()
