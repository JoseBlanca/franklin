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


import logging, os, tempfile, copy, sys
import cPickle as pickle
from itertools import imap, ifilter
from tempfile import gettempdir, NamedTemporaryFile

try:
    import psubprocess
except ImportError:
    pass

import franklin
from franklin.seq.readers import seqs_in_file
from franklin.pipelines.seq_pipeline_steps import SEQPIPELINES, SEQ_STEPS
from franklin.pipelines.snv_pipeline_steps import SNV_PIPELINES, SNV_STEPS
from franklin.pipelines.annotation_steps import ANNOT_STEPS
from franklin.seq.readers import guess_seq_file_format
from franklin.seq.writers import SequenceWriter
from franklin.utils.misc_utils import DisposableFile

# Join the pipelines in PIPELINE
PIPELINES = dict(SEQPIPELINES.items() + SNV_PIPELINES.items())

STEPS = SEQ_STEPS
STEPS.extend(SNV_STEPS)
STEPS.extend(ANNOT_STEPS)

def _get_name_in_config(step):
    'It returns the name in config or the name'
    return step.get('name_in_config', step['name'])

def configure_pipeline(pipeline, configuration):
    '''It chooses the proper pipeline and configures it.'''
    #only for the tests, this function should accept only lists, strs won't be
    #supported
    if isinstance(pipeline, str):
        pipeline = PIPELINES[pipeline]

    # This is done to be able to use the same step more than once. For that we
    # need to have indexed the configuration in different names but knowing wich
    # is the pipeline step

    # set the configuration in the pipeline
    for step in pipeline:
        name_in_config = _get_name_in_config(step)
        if name_in_config in configuration:
            for key, value in configuration[name_in_config].items():
                step['arguments'][key] = value

#    # Here I check that none of the arguments have a none value
#    for step in pipeline:
#        for key, value in step['arguments'].items():
#            # If the step is remove_adaptors, the vectors value can be None
#            if (value is None and
#                (key == 'vectors' and step['name'] != 'remove_adaptors') or :
#
#                msg = 'Parameter %s in step %s from pipeline %s must be set' % \
#                            (key, step['name'], pipeline)
#                raise RuntimeError(msg)
    return pipeline

def _get_func_tools(processes):
    'It returns a mapping function from multiprocessing or itertools'
    if isinstance(processes, bool) and not processes:
        map_ = imap
    else:
        #multiprocessing
        import multiprocessing
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
    if pipeline is None:
        return items
    # We configure the pipeline depending on the sequences type and
    # configuration parameters
    pipeline_steps = configure_pipeline(pipeline, configuration)

    #we create all the cleaner functions
    cleaner_functions = {}
    for analysis_step in pipeline_steps:
        function_factory = analysis_step['function']
        if analysis_step['arguments']:
            arguments = analysis_step['arguments']
        else:
            arguments = None

        msg = "Performing: %s" % analysis_step['comment']
        logging.info(msg)

        # Create function adding parameters if they need them
        if arguments is None:
            cleaner_function = function_factory()
        else:
            #pylint:disable-msg=W0142
            cleaner_function = function_factory(**arguments) #IGNORE:W0142
        cleaner_functions[_get_name_in_config(analysis_step)] = cleaner_function

    #are we multiprocessing?
    functs = _get_func_tools(processes)

    #now use use the cleaner functions using the mapper functions
    for analysis_step in pipeline_steps:
        step_name = _get_name_in_config(analysis_step)
        cleaner_function = cleaner_functions[step_name]
        type_ = analysis_step['type']
        if type_ == 'mapper':
            filtered_items = functs['map'](cleaner_function, items)
        elif type_ == 'filter':
            filtered_items = functs['filter'](cleaner_function, items)
        elif type_ == 'bulk_processor':
            filtered_items = cleaner_function(items)
        items = filtered_items

        msg = "Analysis step prepared: %s" % analysis_step['comment']
        logging.info(msg)

    logging.info('Done!')

    return items

def process_sequences_for_script(in_fpath_seq, file_format,
                                 pipeline, configuration, out_fpath):
    '''It returns a repr file with the processed sequences

    The pipeline and configuration should be pickled object.
    '''
    pipeline = pickle.loads(pipeline)
    configuration = pickle.loads(configuration)
    #the pipeline is now a list of strs we should convert it into a list of
    #dicts
    steps = dict([(step['name'], step) for step in STEPS])
    #the pipeline step do not have real functions because they were not
    #pickeable, so we have to put the functions again in the steps
    for step in pipeline:
        step_name = step['name']
        step['function'] = steps[step_name]['function']

    processed_seqs = _process_sequences(open(in_fpath_seq), in_fhand_qual=None,
                                        file_format=file_format,
                                        pipeline=pipeline,
                                        configuration=configuration)
    #now we write all seq in the file
    out_fhand = open(out_fpath, 'a')
    writer = SequenceWriter(fhand=out_fhand,
                            file_format='pickle')
    for sequence in processed_seqs:
        writer.write(sequence)
    out_fhand.close()

def _parallel_process_sequences(in_fhand_seqs, in_fhand_qual, file_format,
                                pipeline, configuration, processes):
    '''It returns a generator with the processed sequences

    This function does the job calling an external script that does the
    processing. For this calling uses psubprocess doing in fact the
    parallelization
    '''
    if in_fhand_qual:
        #a splitter for both seq and qual should be used
        raise NotImplementedError
    #we have to transform the pipeline list into a list of strs, otherwise
    #it wouldn't be possible to send them to an external script.
    str_pipeline = []
    for step in pipeline:
        str_step = copy.copy(step)
        str_step['function'] = 'function'
        str_pipeline.append(step)

    #everything should be pickle because it will run with an external script
    pipeline = pickle.dumps(pipeline)
    configuration = pickle.dumps(configuration)

    #we need a file that will be removed when the close is called in it
    fhand, out_fpath = tempfile.mkstemp()
    os.close(fhand)

    #debug = 'function'
    #debug = 'subprocess'
    debug = False
    if debug == 'function':
        process_sequences_for_script(in_fhand_seqs.name, file_format,
                                     pipeline, configuration, out_fpath)
    else:
        cmd = [sys.executable]
        cmd2 = os.path.join(franklin.__path__[0], 'process_sequences.py')
        cmd.append(os.path.abspath(cmd2))
        cmd.extend([in_fhand_seqs.name,
                    file_format, pipeline, configuration, out_fpath,
                    gettempdir()])
        if file_format == 'fasta':
            splitter = '>'
        elif file_format in ('fastq', 'sfastq', 'ifastq'):
            splitter = 'fastq'
        elif file_format == 'repr':
            splitter = 'SeqWithQual'
        elif file_format in ('json', 'pickle'):
            splitter = 'blank_line'
        else:
            if file_format:
                msg = 'No parallel splitter for format ' + file_format
            else:
                msg = 'A file with an unknown format cannot be split'
            raise NotImplementedError(msg)
        cmd_def = [{'options': 2, 'io': 'in', 'splitter':splitter},
                   {'options':-2, 'io': 'out'}]
        stdout = NamedTemporaryFile()
        stderr = NamedTemporaryFile()
        if debug == 'subprocess':
            import subprocess
            process = subprocess.Popen(cmd, stdout=stdout,
                                        stderr=stderr)
        else:
            process = psubprocess.Popen(cmd, cmd_def=cmd_def,
                                        stdout=stdout,
                                        stderr=stderr,
                                        splits=processes)
        retcode = process.wait()
        if retcode:
            stdout = open(stdout.name).read()
            stderr = open(stderr.name).read()
            msg = 'Running process seqs script failed.\n stdout: %s\n stderr: %s\n' % (stdout, stderr)
            raise RuntimeError(msg)
        stdout.close()
        stderr.close()
    return seqs_in_file(DisposableFile(out_fpath), format='pickle')

def _process_sequences(in_fhand_seqs, in_fhand_qual, file_format, pipeline,
                                          configuration):
    'It returns a generator with the processed sequences'
    sequences = seqs_in_file(in_fhand_seqs, in_fhand_qual, file_format)

    # the pipeline that will process the generator is build
    processed_seqs = _pipeline_builder(pipeline, sequences, configuration)
    return processed_seqs

def seq_pipeline_runner(pipeline, configuration, in_fhands, file_format=None,
                        writers=None, processes=False):

    '''It runs all the analysis for the given sequence pipeline.

    It takes one or two input files and one or two output files. (Fasta files
    with the sequence and quality).
    A working directory can be given in which the analysis intermediate files
    will be created. If not given a temporary directory will be created that
    will be removed once the analysis is completed.
    If the checkpoints are requested an intermediate file for every step will be
    created.
    '''
    if isinstance(pipeline, str):
        pipeline = PIPELINES[pipeline]

    if file_format is None:
        file_format = guess_seq_file_format(in_fhands['in_seq'])

    # Here we extract our input/output files
    in_fhand_seqs = in_fhands['in_seq']
    if 'in_qual' in in_fhands:
        in_fhand_qual = in_fhands['in_qual']
    else:
        in_fhand_qual = None

    # Here the SeqRecord generator is created
    processes = None if processes == 1 else processes
    if processes:
        sequences = _parallel_process_sequences(in_fhand_seqs,
                                                in_fhand_qual,
                                                file_format, pipeline,
                                                configuration,
                                                processes)
    else:
        sequences = _process_sequences(in_fhand_seqs, in_fhand_qual,
                                       file_format, pipeline,
                                       configuration)

    # The SeqRecord generator is consumed
    for sequence in sequences:
        for writer in writers.values():
            writer.write(sequence)

    # Some of the writers needs to close in order to finish its work
    feature_counter = {}
    for wtype, writer in writers.items():
        if 'close' in dir(writer):
            writer.close()
        feature_counter[wtype] = writer.num_features

    return feature_counter
