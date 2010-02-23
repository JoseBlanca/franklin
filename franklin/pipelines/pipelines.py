'''This module holds some utilities to build sequence cleaning pipelines.

Several pipelines are defined and they can be run using the function
_pipeline_builder. A pipeline consists of a list of several steps that a run
sequentially by the _pipeline_builder. It is very easy to create new pipeline
because they're just a list of dicts that include the name of a function. For
instance there are pipelines defined to clean short and long sequences, to
mask them using repeat masker, etc.

The pipeline can hold three step types: filter, mapper and bulk_processor. They
differ in the interface of the function that process the sequences:
    - mapper: These functions take a sequence and return a new processed
    sequence.
    - filter: These functions take a sequence and return True or False according
    to the match of some criteria by the sequence.
    - bulk_processor: These functions take a sequence iterator and return a new
    sequence iterator will the processed sequence.
The pipeline runner knows how to use these three kinds of steps to filter and
modify the sequences.
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
from itertools import imap, ifilter
import multiprocessing
from franklin.seq.readers import seqs_in_file
from franklin.pipelines.seq_pipeline_steps import SEQPIPELINES
from franklin.pipelines.snv_pipeline_steps import SNV_PIPELINES
from franklin.seq.readers import guess_seq_file_format
from franklin.seq.writers import SequenceWriter
from franklin.snv.writers import VariantCallFormatWriter


# Join the pipelines in PIPELINE
PIPELINES = dict(SEQPIPELINES.items() + SNV_PIPELINES.items())


def configure_pipeline(pipeline, configuration):
    '''It chooses the proper pipeline and configures it.'''

    if isinstance(pipeline, str):
        seq_pipeline  = PIPELINES[pipeline]
    else:
        seq_pipeline = pipeline

    # set the configuration in the pipeline
    for step in seq_pipeline:
        step_name = step['name']
        if step_name in configuration:
            for key, value in configuration[step_name].items():
                step['arguments'][key] = value

    # Here I check that none of the arguments have a none value
    for step in seq_pipeline:
        for key, value in step['arguments'].items():
            if value is None:
                msg = 'Parameter %s in step %s from pipeline %s must be set' % \
                            (key, step['name'], pipeline)
                raise RuntimeError(msg)
    return seq_pipeline

def _get_func_tools(processes):
    'It returns a mapping function from multiprocessing or itertools'
    if isinstance(processes, bool) and not processes:
        map_ = imap
    else:
        #multiprocessing
        if not isinstance(processes, int):
            processes = None    #number determined by cpu_count
        pool = multiprocessing.Pool(processes)
        map_ = pool.imap_unordered
        raise RuntimeError('Multiprocessing not supported yet')
    filter_ = ifilter
    return {'map': map_, 'filter': filter_}

def _pipeline_builder(pipeline, items, configuration=None, processes=False):
    '''It runs all the analysis for the given pipeline.

    It takes one or two input files and one or two output files. (Fasta files
    with the sequence and quality).
    A working directory can be given in which the analysis intermediate files
    will be created. If not given a temporary directory will be created that
    will be removed once the analysis is completed.
    If the checkpoints are requested an intermediate file for every step will be
    created.
    '''
    if configuration is None:
        configuration = {}
    # We configure the pipeline depending on the sequences type and
    # configuration parameters
    pipeline_steps = configure_pipeline(pipeline, configuration)

    # List of temporary files created by the bulk processors.
    # we need to keep them until the analysis is done because some seq_iters
    # may depend on them
    temp_bulk_files = []

    #are we multiprocessing?
    functs = _get_func_tools(processes)

    for analisis_step in pipeline_steps:
        function_factory  = analisis_step['function']
        type_     = analisis_step['type']
        if analisis_step['arguments']:
            arguments = analisis_step['arguments']
        else:
            arguments = None

        msg = "Performing: %s" % analisis_step['comment']
        logging.info(msg)
        #print (msg)

        # Crete function adding parameters if they need them
        if arguments is None:
            cleaner_function = function_factory()
        else:
            #pylint:disable-msg=W0142
            cleaner_function = function_factory(**arguments) #IGNORE:W0142

        if type_ == 'mapper':
            filtered_items = functs['map'](cleaner_function, items)
        elif type_ == 'filter':
            filtered_items   = functs['filter'](cleaner_function, items)
        elif type_ == 'bulk_processor':
            filtered_items, fhand_outs = cleaner_function(items)
            temp_bulk_files.append(fhand_outs)

        items = filtered_items

        #more logging
        msg = "Finished: %s" % analisis_step['comment']
        logging.info(msg)

    temp_bulk_files = None
    logging.info('Done!')

    return items

WRITERS = {'sequence': SequenceWriter,
           'vcf': VariantCallFormatWriter}

def seq_pipeline_runner(pipeline, configuration, io_fhands, file_format=None,
                        processes=False):

    '''It runs all the analysis for the given sequence pipeline.

    It takes one or two input files and one or two output files. (Fasta files
    with the sequence and quality).
    A working directory can be given in which the analysis intermediate files
    will be created. If not given a temporary directory will be created that
    will be removed once the analysis is completed.
    If the checkpoints are requested an intermediate file for every step will be
    created.
    '''
    if file_format is None:
        file_format = guess_seq_file_format(io_fhands['in_seq'])

    # Here we extract our input/output files
    in_fhand_seqs  = io_fhands['in_seq']
    if 'in_qual' in io_fhands:
        in_fhand_qual  = io_fhands['in_qual']
    else:
        in_fhand_qual  = None

    # Here the SeqRecord generator is created
    sequences = seqs_in_file(in_fhand_seqs, in_fhand_qual, file_format)

    # the pipeline that will process the generator is build
    filtered_seq_iter = _pipeline_builder(pipeline, sequences, configuration,
                                        processes)

    #which outputs do we want?
    writers = []
    for output, fhand in io_fhands['outputs'].items():
        if output == 'quality':
            pass
        else:
            writer_klass = WRITERS[output]
            if output == 'sequence':
                if 'quality' in io_fhands['outputs']:
                    qual_fhand = io_fhands['outputs']['quality']
                else:
                    qual_fhand = None
                writer = writer_klass(fhand=fhand,
                                         qual_fhand=qual_fhand,
                                         file_format=file_format)
            elif output == 'repr':
                file_format = 'repr'
                writer = writer_klass(fhand=fhand, file_format='repr')
            elif output == 'vcf':
                ref_name = os.path.basename(io_fhands['in_seq'].name)
                writer = writer_klass(fhand=fhand, reference_name=ref_name)
            else:
                writer = writer_klass(fhand)
            writers.append(writer)

    # The SeqRecord generator is consumed
    for sequence in filtered_seq_iter:
        for writer in writers:
            writer.write(sequence)
