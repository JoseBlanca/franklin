'''
Created on 03/12/2009

@author: peio
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

from franklin.seq.seq_cleaner import (create_vector_striper,
                                      create_adaptor_striper,
                                      create_striper_by_quality,
                                      create_striper_by_quality_lucy,
                                      create_striper_by_quality_trimpoly,
                                      create_masker_for_polia,
                                      create_masker_for_low_complexity,
                                      create_re_word_striper,
                                      create_edge_stripper, create_upper_mapper,
                                      create_seq_trim_and_masker,
                                      create_double_encoding_mapper)

from franklin.seq.seq_filters import (create_length_filter,
                                      create_solid_quality_filter,
                                      create_similar_seqs_filter)

filter_similar_seqs = {'function':create_similar_seqs_filter,
           'arguments':{'db': None, 'blast_program':None},
           'type': 'filter',
           'name': 'filter_similar_seqs',
           'comment': 'It filters similar seqs from a reads database'}

double_encoding = {'function':create_double_encoding_mapper,
           'arguments':{},
           'type': 'mapper',
           'name': 'double_coding',
           'comment': 'It makes double coding to colorspace sequences'}

up_case = {'function':create_upper_mapper,
           'arguments':{},
           'type': 'mapper',
           'name': 'up_case',
           'comment': 'It convers the sequence to upper case'}

#pylint:disable-msg=C0103
remove_vectors_blastdb = {'function':create_vector_striper,
                          'arguments':{'vectors':None,
                                       'vectors_are_blastdb':True},
                          'type': 'mapper',
                          'name': 'remove_vectors_blastdb',
                          'comment': 'Remove vector using vector db'}
remove_vectors_file = {'function':create_vector_striper,
                       'arguments':{'vectors':None,
                                    'vectors_are_blastdb':False},
                       'type': 'mapper',
                       'name': 'remove_vectors_file',
                       'comment': 'Remove vector using vector db'}

remove_adaptors = {'function':create_adaptor_striper,
                   'arguments':{'adaptors':None},
                   'type': 'mapper',
                   'name': 'remove_adaptors',
                   'comment': 'Remove adaptors'}

strip_quality = {'function': create_striper_by_quality,
                      'arguments':{'quality_treshold':20,
                                   'quality_window_width':1,
                                   'only_3_end':False},
#min_quality_bases=None, min_seq_length=None, quality_window_width=None },
                      'type':'mapper',
                      'name':'strip_quality',
                      'comment':'Strip low quality with our algorithm'}
strip_quality_3 = {'function': create_striper_by_quality,
                      'arguments':{'quality_treshold':20,
                                   'quality_window_width':1,
                                   'only_3_end':True},
#min_quality_bases=None, min_seq_length=None, quality_window_width=None },
                      'type':'mapper',
                      'name':'strip_quality',
                    'comment':"Strip low quality from 3 end with our algorithm"}

strip_quality_lucy = {'function': create_striper_by_quality_lucy,
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

filter_short_seqs = {'function': create_length_filter,
                     'arguments':{'length':30, 'count_masked': False},
                     'type':'filter' ,
                     'name':'remove_short',
                     'comment': 'Remove seq shorter than X nt'}
sequence_trimmer = {'function': create_seq_trim_and_masker,
                     'arguments':{},
                     'type':'mapper' ,
                     'name':'trim_and_mask_seq',
                     'comment': 'Trim and mask sequences'}


edge_remover = {'function':create_edge_stripper,
                  'arguments':{},
                  'type': 'mapper',
                  'name': 'edge_removal',
                  'comment': 'Strip given edge lengths. Both sides'}
# words
remove_short_adaptors = {'function': create_re_word_striper,
                         'arguments' : {'words':None},
                         'type'      : 'mapper',
                         'name'      : 'remove_short_adaptors',
                   'comment'   : 'It removes the given regexs from the sequence'
              }

solid_quality = {'function': create_solid_quality_filter,
                'arguments' : {},
                'type'      : 'filter',
                'name'      : 'solid_quality',
                'comment'   : 'It filter by solid quality'
                }


################################################################################
# PIPELINES
################################################################################

SEQPIPELINES = {
    'sanger_with_qual'   : [up_case, remove_adaptors, strip_quality_lucy,
                            remove_vectors_blastdb, remove_vectors_file,
                            remove_adaptors, mask_low_complexity,
                            remove_short_adaptors, edge_remover,
                            sequence_trimmer, filter_short_seqs],

    'sanger_without_qual': [up_case, remove_vectors_blastdb,
                            remove_vectors_file, remove_adaptors,
                            strip_quality_by_n,
                            mask_low_complexity, remove_short_adaptors,
                            edge_remover, sequence_trimmer, filter_short_seqs],

    'solexa'             : [up_case, remove_adaptors, strip_quality,
                            sequence_trimmer, filter_short_seqs],

    'adaptors'           : [remove_adaptors, sequence_trimmer,
                            filter_short_seqs],

    'mask_dust'          : [mask_polia, mask_low_complexity, sequence_trimmer],

    'word_masker'        : [remove_short_adaptors, sequence_trimmer,
                            filter_short_seqs],

    'solid'              : [double_encoding, solid_quality, strip_quality_3,
                            sequence_trimmer, filter_short_seqs]}

SEQ_STEPS = [remove_vectors_blastdb, remove_vectors_file, remove_adaptors,
             strip_quality, strip_quality_lucy, strip_quality_by_n,
             strip_quality_by_n, mask_polia, mask_low_complexity,
             sequence_trimmer, filter_short_seqs, edge_remover,
             remove_short_adaptors, up_case, solid_quality, strip_quality_3,
             filter_similar_seqs]
