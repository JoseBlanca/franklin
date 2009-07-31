'''This module holds some utilities to build sequence cleaning pipelines.

Several pipelines are defined and they can be run using the function
pipeline_runner. A pipeline consists of a list of several steps that a run
sequentially by the pipeline_runner. It is very easy to create new pipeline
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


import os, biolib, logging
from itertools import imap, ifilter

from biolib.seq_cleaner import (create_vector_striper_by_alignment,
                                create_striper_by_quality,
                                create_striper_by_quality_lucy,
                                create_striper_by_quality_lucy2,
                                create_striper_by_quality_trimpoly,
                                create_masker_for_polia,
                                create_masker_for_low_complexity,
                                create_masker_repeats_by_repeatmasker)
from biolib.contig_cleaner import (create_contig_read_stripper,
                                   create_read_number_contig_filter,
                                   create_non_matched_region_stripper)
from biolib.snp_cleaner import (create_bad_quality_allele_remover,
                                create_cap_enzyme_filter,
                                create_pic_filter,
                                create_second_allele_number_filter,
                                create_seqvar_close_to_limit_filter)

from biolib.seq_filters        import create_length_filter
from biolib.biolib_seqio_utils import seqs_in_file, write_fasta_file


DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

################################################################################
# PIPELINE CLEANING STEPS
################################################################################

####  Seqs cleaning #################

#pylint:disable-msg=C0103
remove_vectors = {'function':create_vector_striper_by_alignment,
                  'arguments':{'vectors':None, 'aligner':'blast'},
                  'type': 'mapper',
                  'name': 'remove_vectors',
                  'comment': 'Remove vector using vector db'}



remove_adaptors_solexa = {'function':create_vector_striper_by_alignment,
       'arguments':{'vectors':os.path.join(DATA_DIR, 'standar_solexa_adaptors'),
                     'aligner':'exonerate'},
       'type': 'mapper',
       'name': 'remove_adaptors',
       'comment': 'Remove our adaptors'}

strip_quality = {'function': create_striper_by_quality,
                      'arguments':{'quality_treshold':20,
                                   'quality_window_width':1},
#min_quality_bases=None, min_seq_length=None, quality_window_width=None },
                      'type':'mapper',
                      'name':'strip_quality',
                      'comment':'Strip low quality with our algorithm'}

strip_quality_lucy = {'function': create_striper_by_quality_lucy,
                      'arguments':{},
                      'type':'mapper',
                      'name':'strip_lucy',
                      'comment':'Strip low quality with lucy'}

strip_quality_lucy2 = {'function': create_striper_by_quality_lucy2,
                      'arguments':{},
                      'type':'bulk_processor',
                      'name':'strip_lucy',
                      'comment':'Strip low quality with lucy'}

strip_quality_by_n = {'function': create_striper_by_quality_trimpoly,
                          'arguments': {},
                          'type':'mapper',
                          'name':'strip_trimpoly',
                          'comment':'Strip low quality with trimpoly'}

mask_polia         = {'function': create_masker_for_polia,
                       'arguments': {},
                       'type':'mapper',
                       'name':'mask_polia',
                       'comment':'Mask poli A regions'}


mask_low_complexity = {'function': create_masker_for_low_complexity,
                       'arguments': {},
                       'type':'mapper',
                       'name':'mask_low_complex',
                       'comment':'Mask low complexity regions'}

mask_repeats = {'function':create_masker_repeats_by_repeatmasker ,
                'arguments':{'species':'eudicotyledons'},
                'type': 'mapper',
                'name': 'mask_repeats',
                'comment':'Mask repeats with repeatmasker'}

filter_short_seqs_sanger = {'function': create_length_filter,
                     'arguments':{'length':100, 'count_masked': False},
                     'type':'filter' ,
                     'name':'remove_short',
                     'comment': 'Remove seq shorter than 100 nt'}

filter_short_seqs_solexa = {'function': create_length_filter,
                            'arguments':{'length':22, 'count_masked': False},
                            'type':'filter' ,
                            'name':'remove_short',
                            'comment': 'Remove seq shorter than 22 nt'}

## Contig cleaning   ####

contig_extreme_stripper = {'function':create_contig_read_stripper,
                        'arguments':{'length_to_strip':12},
                        'type': 'mapper',
                        'name': 'Contig_stripper',
                        'comment':'It strips reads extremes in contig extremes'}

contig_non_matched_stripper = {'function':create_non_matched_region_stripper,
                        'arguments':{},
                        'type': 'mapper',
                        'name': '',
                        'comment':''}

contig_read_num_filter = {'function':create_read_number_contig_filter,
                        'arguments':{'min_read_number':4 },
                        'type': 'filter',
                        'name': 'filter_by_read_number',
                        'comment':'It filters by read number in contig'}


# Snp cleaning /filtering ####

snp_remove_baq_quality_alleles = {'function':create_bad_quality_allele_remover,
                                  'arguments':{'qual_threshold':20},
                                  'type':'mapper',
                                  'name':'bad_quality_allele_striper',
                                  'comment': 'It removes bad quality alleles'}
snp_second_allele_filter = {'function':create_second_allele_number_filter,
                            'arguments':{'number_2allele':2},
                            'type':'filter',
                            'name':'second_allele_num',
                            'comment': 'It filters by second allele number'}

snp_limit_filter = {'function':create_seqvar_close_to_limit_filter,
                    'arguments':{'max_distance':12},
                    'type':'filter',
                    'name':'limit_distance',
                    'comment': 'It filters by the distancie to the limit'}

pic_filter = {'function':create_pic_filter,
              'arguments':{'min_pic': 0.05},
              'type':'filter',
              'name':'pic_filter',
              'comment': 'It filters the snp by its pic calcule'}

cap_enzyme_filter  = {'function':  create_cap_enzyme_filter,
                      'arguments': {'all_enzymes':True},
                      'type':      'filter',
                      'name':      'enzyme_filter',
                      'comment':  'It filters by enzyme'}



################################################################################
# PIPELINES
################################################################################

PIPELINES = {'sanger_with_qual' : [remove_vectors, strip_quality_lucy2,
                                mask_low_complexity, filter_short_seqs_sanger ],

            'sanger_without_qual': [remove_vectors, strip_quality_by_n,
                               mask_low_complexity, filter_short_seqs_sanger ],

            'repeatmasker' : [mask_repeats, filter_short_seqs_sanger],

            'solexa'       : [remove_adaptors_solexa, mask_low_complexity,
                           mask_polia, strip_quality, filter_short_seqs_solexa],
            'contig_clean' : [contig_extreme_stripper,
                              contig_non_matched_stripper,
                              contig_read_num_filter],
         'snp_clean':[snp_remove_baq_quality_alleles, snp_second_allele_filter,
                      snp_limit_filter], #cap_enzyme_filter, pic_filter]
            }
################################################################################

def configure_pipeline(pipeline, configuration):
    '''It chooses the proper pipeline and configures it.'''

    seq_pipeline  = PIPELINES[pipeline]

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

def pipeline_runner(pipeline, items, configuration=None):
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
    # configuratiom parameters
    pipeline_steps = configure_pipeline(pipeline, configuration)

    # List of temporary files created by the bulk processors.
    # we need to keep them until the analysis is done because some seq_iters
    # may depend on them
    temp_bulk_files = []

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
            # pylint:disable-msg=W0142
            cleaner_function = function_factory(**arguments)

        if type_ == 'mapper':
            filtered_items = imap(cleaner_function, items)
        elif type_ == 'filter':
            filtered_items   = ifilter(cleaner_function, items)
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


def seq_pipeline_runner(pipeline, configuration, io_fhands, file_format=None):

    '''It runs all the analysis for the given sequence pipeline.

    It takes one or two input files and one or two output files. (Fasta files
    with the sequence and quality).
    A working directory can be given in which the analysis intermediate files
    will be created. If not given a temporary directory will be created that
    will be removed once the analysis is completed.
    If the checkpoints are requested an intermediate file for every step will be
    created.
    '''
    # Here we extract our input/output files
    in_fhand_seqs  = io_fhands['in_seq']
    in_fhand_qual  = io_fhands['in_qual']
    out_fhand_seq  = io_fhands['out_seq']
    out_fhand_qual = io_fhands['out_qual']

    # Here starts the analisis
    seq_iter = seqs_in_file(in_fhand_seqs, in_fhand_qual, file_format)

    #run the pipeline
    filtered_seq_iter = pipeline_runner(pipeline, seq_iter, configuration)


    # write result
    write_fasta_file(filtered_seq_iter, out_fhand_seq,  out_fhand_qual)











