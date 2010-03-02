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
from franklin.snv.snv_annotation import (calculate_maf_frequency,
                                         snvs_in_window, calculate_snv_kind,
                                         calculate_cap_enzymes)

FILTER_DESCRIPTIONS = {
    'uniq_contiguous':
        {'id': 'UCR',
         'description':'A blast in the near region gave several matches'},
    'close_to_intron':
        {'id': 'I%2d',
         'description':'An intron is located closer than %2d base pairs'},
    'High_variable_region':
        {'id': 'HVR%2d',
    'description':'The snv is in a region with more than %2d % of variability'},
    }

def _add_filter_result(snv, filter_name, result, threshold=None):
    'It adds the filter to the SeqFeature qualifiers'
    qualifiers = snv.qualifiers
    if 'filters' not in qualifiers:
        qualifiers['filters'] = {}
    if filter_name not in qualifiers['filters']:
        qualifiers['filters'][filter_name] = {}
    qualifiers['filters'][filter_name][threshold] = result

def _get_filter_result(snv, filter_name, threshold=None):
    'It gets the result of a filter. Returns None if the filter is not done'
    qualifiers = snv.qualifiers
    if 'filters' not in qualifiers:
        return None
    try:
        result = qualifiers['filters'][filter_name][threshold]
        return result
    except KeyError:
        return None

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
            # Check if it is already done
            previous_result = _get_filter_result(snv, 'uniq_contiguous',
                                                 threshold=distance)
            if previous_result is not None:
                continue

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

            _add_filter_result(snv, 'uniq_contiguous', result, distance)

    return unique_contiguous_region_filter

def create_close_to_intron_filter(distance):
    '''It returns a filter that filters snv by the proximity to introns.

    The introns should have been annotated before.'''

    def close_to_intron_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'close_to_intron',
                                                 threshold=distance)
            if previous_result is not None:
                continue

            location = snv.location.start.position
            result = True
            for intron in sequence.get_features(kind='intron'):
                if abs(location - intron.location.start.position) < distance:
                    result = False
            _add_filter_result(snv, 'close_to_intron', result,
                               threshold=distance)
    return close_to_intron_filter

def create_high_variable_region_filter(max_variability, window=None):
    'It returns a filter that filters snvs by region variability.'

    def high_variable_region_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        snvs = list(sequence.get_features(kind='snv'))
        for snv in snvs:
            threshold = (max_variability, window)
            previous_result = _get_filter_result(snv, 'high_variable_reg',
                                                 threshold=threshold)
            if previous_result is not None:
                continue
            if window is None:
                snv_num      = len(snvs)
                total_length = len(sequence)
            else:
                total_length = window
                snv_num = snvs_in_window(snv, snvs, window)
            variability = (snv_num/float(total_length)) * 100
            if variability > max_variability:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'high_variable_reg', result,
                               threshold=threshold)
    return high_variable_region_filter




def create_close_to_snv_filter(proximity):
    '''It returns a filter that filters snv by the proximity to other snvs.

    If the snv has another snv closer than DISTANCE, then this snv is
    filtered out'''
    def close_to_snv_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        snvs = list(sequence.get_features(kind='snv'))
        for snv in snvs:
            previous_result = _get_filter_result(snv, 'close_to_snv',
                                                 threshold=proximity)
            if previous_result is not None:
                continue

            num_snvs = snvs_in_window(snv, snvs, proximity * 2)
            if num_snvs > 1:
                result = True
            else:
                result = False

            _add_filter_result(snv, 'close_to_snv', result,
                               threshold=proximity)
    return close_to_snv_filter

def create_snv_close_to_limit_filter(distance):
    '''This function is a function factory. This function is a filter than
     return true if those snvs are or not closer to the limit than max_distance
     '''
    def snv_close_to_reference_limit(sequence):
        '''True if the snv variation is close to the contig limits.

        It checks if the snv is not within the range covered by the
        consensus and if it's close to one of its limits. In both cases it will
        return True.
        '''
        if sequence is None:
            return None
        snvs = list(sequence.get_features(kind='snv'))
        for snv in snvs:
            previous_result = _get_filter_result(snv, 'close_to_limit',
                                                 threshold=distance)
            if previous_result is not None:
                continue
            location = int(str(snv.location.start))
            if location < distance or location + distance > len(sequence):
                result = True
            else:
                result = False
            _add_filter_result(snv, 'close_to_limit', result,
                               threshold=distance)
    return  snv_close_to_reference_limit


def create_major_allele_freq_filter(frequency):
    'It filters the snv in a seq by the frecuency of the more frecuent allele'
    def major_allele_freq_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'close_to_limit',
                                                 threshold=frequency)
            if previous_result is not None:
                continue
            maf = calculate_maf_frequency(snv)
            if maf > frequency:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'maf', result,
                               threshold=frequency)

    return major_allele_freq_filter



def create_kind_filter(kind):
    'It filters the snv by its kind'
    def kind_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'by_kind', threshold=kind)
            if previous_result is not None:
                continue

            kind_ = calculate_snv_kind(snv)
            if kind  == kind_:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'by_kind', result, threshold=kind)
    return kind_filter

def create_cap_enzyme_filter(all_enzymes):
    'It filters the snv looking if it is detectable by rstriction enzymes'
    def cap_enzyme_filter(sequence):
        'The filter'
        if sequence is None:
            return None
        for snv in sequence.get_features(kind='snv'):
            previous_result = _get_filter_result(snv, 'cap_enzymes',
                                                 threshold=all_enzymes)
            if previous_result is not None:
                continue
            enzymes = calculate_cap_enzymes(snv, sequence,
                                            all_enzymes=all_enzymes)
            if len(enzymes) != 0:
                result = True
            else:
                result = False
            _add_filter_result(snv, 'cap_enzymes', result, threshold=all_enzymes)
    return cap_enzyme_filter



