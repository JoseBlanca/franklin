#!/usr/bin/env python
from franklin.seq.seq_cleaner import check_sequences_length

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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

'It cleans reads using the given configration'

import sys
from optparse import OptionParser

from franklin.utils.misc_utils import get_num_threads
from franklin.seq.readers import guess_seq_file_format
from franklin.seq.writers import SequenceWriter
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.pipelines.seq_pipeline_steps import (double_coding,
                                                   sequence_trimmer,
                                                   edge_remover,
                                                   solid_quality,
                                                   strip_quality, up_case,
                                                   remove_adaptors,
                                                   remove_vectors_blastdb,
                                                   remove_vectors_file,
                                                   strip_quality_lucy,
                                                   strip_quality_by_n,
                                                   mask_low_complexity,
                                                   remove_short_adaptors)

class UnknownFormatError(Exception):
    'Special error for guess_seq_file error'
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PsubprocessNotFounError(Exception):
    'Error to use when psubprocess is not found'

class PlatformError(Exception):
    "Error to use when platform and options don't fit"

class AdaptorError(Exception):
    "Error to use when adaptors length is less than min"
    pass


DEFAULT_VALUES = {'lucy_bracket':[10, 0.02],
                  'lucy_window':[50, 0.08, 10, 0.3],
                  'lucy_error':[0.015, 0.015],
                  'qual_threshold':20,
                  'qual_window':1,
                  'solid_qual_length':10,
                  "solid_qual_threshold":15,
                  'solid_disable_missing_call':True}

def parse_options():
    'It parses the command line arguments'
    usage = "usage: %prog -i seq_file -o out_seq_file -p platform [options]"
    parser = OptionParser(usage=usage)
    parser.add_option('-i', '--seq_in', dest='seq_in',
                    help='Input reads sequence file. [Required]')
    parser.add_option('-q', '--qual_in', dest='qual_in', default=None,
                      help='Input reads quality file. [Optional]')
    parser.add_option('-f', '--format', dest="format", default=None,
                      help='Input read file format. [Optional]')
    parser.add_option('-o', '--seq_out', dest='seq_out',
                    help='Output sequence file, [Required]')
    parser.add_option('-u', '--qual_out', dest='qual_out', default=None,
                      help='Output quality file. [Optional]')
    parser.add_option('-g', '--output_format', dest="out_format", default=None,
                      help='Output read file format. [Optional]')
    parser.add_option('-t', '--threads', dest="threads", default=1, type='int',
                      help='Num threads. Choosing 0 it will use all cores.')
    pl_msg = 'Sequencing platform. (sanger, 454, illumina, solid) [Required]'
    parser.add_option('-p', '--platform', dest="platform", help=pl_msg)
    re_msg = 'Words separated by commas to remove from sequences'
    parser.add_option('-r', '--re_word', dest="re_words", default=None,
                      help=re_msg)
    adapt_msg = 'File with adaptors to remove from reads'
    parser.add_option('-a', '--adaptors_file', dest="adaptors_file",
                      default=None, help=adapt_msg)
    parser.add_option('-v', '--vector_file', dest="vector_file", default=None,
                      help='File with vectors to remove from reads')
    parser.add_option('-d', '--vector_db', dest="vector_db", default=None,
                      help='Blast database to use to remove vectors')
    parser.add_option('-m', '--min_len', dest="min_len",
                      help='Minimun length of reads to filter out')
    parser.add_option('-e', '--edge_trim', dest="edge_trim", default=None,
                      help='Nucleotides to trim in the edges')
    parser.add_option('-n', '--n_percent', dest="n_percent", default=None,
                      help='Number of nucleotides to trim in the edges')
    parser.add_option('-l', '--lucy_splice_file', dest="lucy_splice",
                      default=None, help='Lucy splice file')
    parser.add_option('--lucy_bracket', dest="lucy_bracket",
                      default=DEFAULT_VALUES['lucy_bracket'],
                      help='Lucy bracket paramters')
    parser.add_option('--lucy_window', dest="lucy_window",
                      default=DEFAULT_VALUES['lucy_windows'],
                      help='Lucy bracket parameters')
    parser.add_option('--lucy_error', dest="lucy_error",
                      default=DEFAULT_VALUES['lucy_error'], help='Lucy error parameters')
    parser.add_option('--qual_threshold', dest="qual_threshold",
                      default=DEFAULT_VALUES['qual_treshold'], help='Quality threshold parameters')
    parser.add_option('--qual_window', dest="qual_window",
                      default=DEFAULT_VALUES['qual_window'], help='Quality window parameter')
    parser.add_option('--only_3_end', dest="only_3_end", action="store_true",
                      default=None, help='Quality remove only 3 end')
    parser.add_option('--double_encoding', dest="double_encoding",
                      action='store_true', default=False,
                      help='Double encoding for solid')
    parser.add_option('--solid_qual_length', dest="solid_qual_length",
                      default=DEFAULT_VALUES['solid_qual_length'],
                      help="Number of 5' colors to consider to quality filtering")
    parser.add_option('--solid_qual_threshold', dest="solid_qual_threshold",
                      default=DEFAULT_VALUES['solid_qual_threshold'],
                      help="Minimum mean quality allowable for solid reads.")
    parser.add_option('--solid_disable_missing_call',
                      dest="solid_disable_missing_call", action='store_true',
                      default=DEFAULT_VALUES['solid_disable_missing_call'],
                      help="Minimum mean quality allowable for solid reads.")
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.seq_in is None:
        parser.error("Input seqfile is required")
    else:
        seq_in_fhand = open(options.seq_in)

    if options.qual_in is None:
        qual_in_fhand = None
    else:
        qual_in_fhand = open(options.qual_in)
    format_in = options.format

    if options.seq_out is None:
        parser.error("Output seqfile is required")
    else:
        seq_out_fhand = open(options.seq_out, 'w')

    if options.qual_out is None:
        qual_out_fhand = None
    else:
        qual_out_fhand = open(options.qual_out, 'w')
    format_out = options.out_format
    io = {'seq_in_fhand':seq_in_fhand, 'qual_in_fhand':qual_in_fhand,
          'format_in':format_in, 'format_out':format_out,
          'seq_out_fhand':seq_out_fhand, 'qual_out_fhand':qual_out_fhand}
    io['double_encoding']  = options.double_encoding

    threads = options.threads
    platform = options.platform
    parameters = {}
    parameters['platform']         = platform
    parameters['re_word']          = options.re_word
    parameters['adaptors_file']    = options.adaptors_file
    parameters['vector_file']      = options.vector_file
    parameters['vector_db']        = options.vector_db
    parameters['min_len']          = options.min_len
    parameters['edge_trim']        = options.edge_trim
    parameters['n_percent']        = options.n_percent
    parameters['lucy_splice_file'] = options.lucy_splice_file
    parameters['lucy_bracket']     = options.lucy_bracket
    parameters['lucy_window']      = options.lucy_window
    parameters['lucy_error']       = options.lucy_error
    parameters['qual_threshold']   = options.qual_threshold
    parameters['qual_min_length']  = options.qual_min_length
    parameters['qual_window']      = options.qual_window
    parameters['only_3_end']       = options.only_3_end

    return io, parameters, threads

def _prepare_inputs(io):
    'It prepares input files'
    in_seq_fhand  = io['seq_in_fhand']
    in_qual_fhand = io['qual_in_fhand']

    if io['format_in'] is None:
        try:
            io['format_in'] = guess_seq_file_format(in_seq_fhand)
        except ValueError:
            raise UnknownFormatError
    format_in = io['format_in']

    if format_in == 'fasta' and io['qual_in_fhand'] is None:
        seq_with_qual = False
    else:
        seq_with_qual = True

    return ({'in_seq':in_seq_fhand, 'in_qual':in_qual_fhand}, format_in,
            seq_with_qual)

def _prepare_outputs(io):
    '''It prepares the outputs. It returns a writer instance to use with
    pipeline runner. take into account the format'''
    if io['format_out'] is None:
        if io['format_in'] == 'csfasta':
            io['format_out'] = 'sfastq'
        else:
            io['format_out'] = io['format_in']

    writer = SequenceWriter(io['seq_out_fhand'], io['format_out'],
                            qual_fhand=io['qual_out_fhand'])
    return writer

def _get_num_threads(threads):
    'it calculates the threads'
    if threads == 0:
        try:
            import psubprocess
        except ImportError:
            raise PsubprocessNotFounError

        threads = True
    return get_num_threads(threads)

def _prepare_pipeline(parameters, qual):
    'It prepares_the pipeline'
    if parameters['adaptor_file']:
        try:
            check_sequences_length(open(parameters['adaptor_file']), 15)
        except ValueError as msg:
            raise AdaptorError('Adaptor length error:\n%s' % msg)

    platform = parameters['platform']
    if platform == 'solid':
        pipeline, configuration = _prepare_pipeline_for_solid(parameters)
    elif platform == 'illumina':
        pipeline, configuration = _prepare_pipeline_for_illumina(parameters)
    elif platform in ('454', 'sanger') :
        pipeline, configuration = _prepare_pipeline_for_sanger(parameters,
                                                               qual)

    # triming and seq_lengh steps
    pipeline.append(sequence_trimmer)
    if parameters['min_len']:
        pipeline.append(sequence_trimmer)
        configuration['length'] = parameters['min_length']
    return pipeline, configuration

def _prepare_pipeline_for_solid(parameters):
    'it prepares the solid pipeline'
    pipeline      = []
    configuration = {}

    if parameters['double_coding']:
        pipeline.append('double_coding')
        configuration[double_coding] = {}

    #solid quality
    pipeline.append(solid_quality)
    configuration['solid_quality'] = {'length':parameters['solid_qual_length'],
                                      'threshold':parameters['solid_qual_threshold'],
                                      'call_missing':parameters['solid_disable_missing_call']}

    #strip_quality
    pipeline.append(strip_quality)
    configuration['strip_quality'] = {'quality_treshold':parameters['qual_threshold'],
                                      'quality_window_width':parameters['qual_window'],
                                      'only_3_end':True}

    return pipeline, configuration

def _prepare_pipeline_for_illumina(parameters):
    'it prepares the illumina pipeline'
    pipeline      = []
    configuration = {}

    #upcase
    pipeline.append(up_case)
    configuration['up_case'] = {}

    #remove_adaptors
    if parameters['adaptors_file']:
        pipeline.append(remove_adaptors)
        configuration['remove_adaptor'] = {'adaptors':parameters['adaptors_file']}
    # strip_quality
    pipeline.append(strip_quality)
    configuration['strip_quality'] = {'quality_treshold':parameters['qual_threshold'],
                                      'quality_window_width':parameters['qual_window'],
                                      'only_3_end':parameters['only_3_end']}

    return pipeline, configuration

def _prepare_pipeline_for_sanger(parameters, qual):
    'it prepares the illumina pipeline'
    pipeline      = []
    configuration = {}

    #upcase
    pipeline.append(up_case)
    configuration['up_case'] = {}

    #remove_adaptors
    if parameters['adaptors_file']:
        pipeline.append(remove_adaptors)
        configuration['remove_adaptor'] = {'adaptors':parameters['adaptors_file']}

    #vector_db
    if parameters['vector_db']:
        pipeline.append(remove_vectors_blastdb)
        configuration['remove_vectors_blastdb'] = {'vectors':parameters['vector_db']}

    #vector_file
    if parameters['vector_file']:
        pipeline.append(remove_vectors_file)
        configuration['remove_vectors_file'] = {'vectors':parameters['vector_file']}

    #lucy
    if qual:
        pipeline.append(strip_quality_lucy)
        configuration['strip_lucy'] = {'bracket':parameters['lucy_bracket'],
                                       'window':parameters['lucy_window'],
                                       'error':parameters['lucy_error']}
        if parameters['vector_file'] and parameters['lucy_splice_file']:
            configuration['strip_lucy']['vector'] = [parameters['vector_file'],
                                                     parameters['lucy_splice_file']]

    # n_percent
    if not qual:
        pipeline.append(strip_quality_by_n)
        configuration['strip_trimpoly'] = {'n_trim_above_percent':parameters['n_percent']}

    #mask_low_complexity
    pipeline.append(mask_low_complexity)
    configuration['mask_low_complex'] = {}

    #remove_short_adaptors
    if parameters['re_words']:
        pipeline.append(remove_short_adaptors)
        configuration['remove_short_adaptors'] = {'words':parameters['re_words']}

    #edge_remover
    if parameters['edge_trimp']:
        left  = None if parameters['edge_trimp'][0] == 0 else parameters['edge_trimp'][1]
        rigth = None if parameters['edge_trimp'][0] == 0 else parameters['edge_trimp'][1]
        pipeline.append(edge_remover)
        configuration['ege_removal'] = {'left_length':left, 'right_length':rigth}

    return pipeline, configuration

def main():
    'The main part'
    try:
        io, parameters, threads = set_parameters()
    except PlatformError as msg:
        sys.stderr.write(msg)
        sys.exit()

    # input
    try:
        in_fhands, informat, seq_with_qual = _prepare_inputs(io)
    except UnknownFormatError:
        msg = "Can not guess infut file format. Please use -f option"
        sys.stderr.write(msg)
        sys.exit()

    # output
    writer = _prepare_outputs(io)

    # threads
    try:
        threads = _get_num_threads(threads)
    except PsubprocessNotFounError:
        sys.stderr.write('In order to use multiprocessin you need to install psubprocess')
        sys.exit()

    #do we have seqs with quality


    # parameters
    try:
        pipeline, configuration = _prepare_pipeline(parameters, seq_with_qual)
    except AdaptorError as msg:
        sys.stderr.write(msg)
        sys.exit()


    # the runner
    seq_pipeline_runner(pipeline, configuration, in_fhands,
                        file_format=informat, writers=[writer],
                        processes=threads)

if __name__ == '__main__':
    main()
