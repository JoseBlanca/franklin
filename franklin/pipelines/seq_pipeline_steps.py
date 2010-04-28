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

from franklin.seq.seq_cleaner import (create_vector_striper_by_alignment,
                                    create_striper_by_quality,
                                    create_striper_by_quality_lucy,
                                    create_striper_by_quality_lucy2,
                                    create_striper_by_quality_trimpoly,
                                    create_masker_for_polia,
                                    create_masker_for_low_complexity,
                                    create_masker_repeats_by_repeatmasker,
                                    create_word_remover,
                                    create_edge_stripper)

from franklin.seq.seq_filters import create_length_filter

#pylint:disable-msg=C0103
remove_vectors = {'function':create_vector_striper_by_alignment,
                  'arguments':{'vectors':None, 'aligner':'blast'},
                  'type': 'mapper',
                  'name': 'remove_vectors',
                  'comment': 'Remove vector using vector db'}

remove_adaptors = {'function':create_vector_striper_by_alignment,
       'arguments':{'vectors':None, 'aligner':'exonerate'},
                    #os.path.join(DATA_DIR, 'standar_solexa_adaptors.fasta'),
       'type': 'mapper',
       'name': 'remove_adaptors',
       'comment': 'Remove adaptors'}

strip_quality = {'function': create_striper_by_quality,
                      'arguments':{'quality_treshold':20,
                                   'quality_window_width':1},
#min_quality_bases=None, min_seq_length=None, quality_window_width=None },
                      'type':'mapper',
                      'name':'strip_quality',
                      'comment':'Strip low quality with our algorithm'}

strip_quality_lucy = {'function': create_striper_by_quality_lucy,
                      'arguments':{'vector':None,'splice_site':None },
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

edge_remover = {'function':create_edge_stripper,
                  'arguments':{},
                  'type': 'mapper',
                  'name': 'edge_removal',
                  'comment': 'Strip given edge lengths. Both sides'}
# words
remove_words = {'function'  : create_word_remover,
              'arguments' : {'words':None},
              'type'      : 'mapper',
              'name'      : 'word_remover',
      'comment'   : 'It removes  the given words in the beginning of a sequence'
              }


################################################################################
# PIPELINES
################################################################################

SEQPIPELINES = {
    'sanger_with_qual'   : [remove_adaptors, strip_quality_lucy2,
                            remove_vectors, mask_low_complexity,
                            remove_words, edge_remover,
                            filter_short_seqs_sanger],

    'sanger_without_qual': [remove_vectors, strip_quality_by_n,
                            mask_low_complexity, remove_words, edge_remover,
                            filter_short_seqs_sanger],

    'repeatmasker'       : [mask_repeats, filter_short_seqs_sanger],

    'solexa'             : [remove_adaptors, strip_quality,
                            filter_short_seqs_solexa],

    'adaptors'           : [remove_adaptors, filter_short_seqs_sanger],

    'mask_dust'          : [mask_polia, mask_low_complexity],

    'word_masker'        : [remove_words, filter_short_seqs_solexa]}


