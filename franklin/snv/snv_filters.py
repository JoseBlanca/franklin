'''
Created on 19/02/2010

@author: jose
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

from Bio import SeqIO

from franklin.utils.cmd_utils import create_runner
from franklin.alignment_search_result import (FilteredAlignmentResults,
                                              get_alignment_parser)
from franklin.seq.seq_analysis import (infer_introns_for_cdna,
                                       similar_sequences_for_blast)
from franklin.seq.readers import guess_seq_file_format

FILTER_DESCRIPTIONS = {
    'uniq_contiguous':
        {'id': 'UCR',
         'description':'A blast in the near region gave several matches'},
    'close_to_intron':
        {'id': 'I%2d',
         'description':'An intron is located closer than %2d base pairs'},
         }

def _add_filter_result(snv, filter_name, result, threshold=None):
    'It adds the filter to the SeqFeature qualifiers'
    qualifiers = snv.qualifiers
    if 'filters' not in qualifiers:
        qualifiers['filters'] = {}
    res = {'result': result}
    if threshold is not None:
        res['threshold'] = threshold
    qualifiers['filters'][filter_name] = res

def create_unique_contiguous_region_filter(distance, genomic_db,
                                           genomic_seqs_fhand):
    '''It returns a filter that removes snv in a region that give more than one
    match or more than one match_parts'''
    parameters = {'database': genomic_db, 'program':'blastn'}
    blast_runner = create_runner(tool='blast', parameters=parameters)
    blast_parser = get_alignment_parser('blast')
    filters = [{'kind'           : 'min_scores',
                'score_key'      : 'similarity',
                'min_score_value': 90,
               },
               {'kind'           : 'min_length',
                'min_length_bp'  : 20,
               }
              ]
    genomic_seqs_index = SeqIO.index(genomic_seqs_fhand.name,
                                     guess_seq_file_format(genomic_seqs_fhand))

    def unique_contiguous_region_filter(sequence):
        '''It filters out the snv in regions repeated in the genome or
        discontiguous'''
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            #we make a blast
            #with the sequence around the snv
            location = snv.location.start.position
            start = location - distance
            end   = location + distance
            if start < 0:
                start = 0
            #print start, end
            seq_fragment = sequence[start:end]
            blast_fhand = blast_runner(seq_fragment)['blast']
            #now we parse the blast
            blast_result = blast_parser(blast_fhand)
            alignments = FilteredAlignmentResults(match_filters=filters,
                                              results=blast_result)
            #are there any similar sequences?
            try:
                alignment = alignments.next()
            except StopIteration:
                #if there is no similar sequence we assume that is unique
                result = True
            #how many matches, it should be only one
            num_hits = len(alignment['matches'])

            if num_hits > 1:
                result = False
            else:
                #how many match parts have the first match?
                #we could do it with the blast result, but blast is not very
                #good aligning, so we realign with est2genome
                blast_fhand.seek(0)
                sim_seqs = similar_sequences_for_blast(blast_fhand)
                introns = infer_introns_for_cdna(sequence=seq_fragment,
                                          genomic_seqs_index=genomic_seqs_index,
                                          similar_sequence=sim_seqs[0])
                if not introns:
                    result = True
                else:
                    result = False

            _add_filter_result(snv, 'uniq_contiguous', result)

    return unique_contiguous_region_filter

def create_close_to_intron_filter(distance):
    '''It returns a filter that filters snv by the proximity to introns.

    The introns should have been annotated before.'''

    def close_to_intron_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            location = snv.location.start.position
            result = True
            for intron in sequence.get_features(kind='intron'):
                if abs(location - intron.location.start.position) < distance:
                    result = False
            _add_filter_result(snv, 'close_to_intron', result,
                               threshold=distance)
    return close_to_intron_filter
