'''This module holds utilities to filter sequences.

The filtering can be divided in two kinds of functions, the ones that return
a bool meaning if the sequence should be removed or not and the ones that
return a masked sequence or trimmed sequence. In this latter case the filter
can also return None if no sequence is left after the filtering process.
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

from biolib.biolib_cmd_utils import  create_runner
from biolib.biolib_utils import temp_multi_fasta_file
from biolib.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)



def create_aligner_filter(aligner_cmd, cmd_parameters, match_filters=None,
                            result_filters=None):
    '''A function factory factory that creates aligner filters.

    It returns a function that will accept a sequence and it will return
    True or False depending on the exonerate outcome.
    parameters is a dictionary and key are defined in ExonerateRunner.
    Required is only the target fasta file
    '''
    #runners = {'blast':BlastRunner, 'exonerate':ExonerateRunner}

    parser = get_alignment_parser(aligner_cmd)
    binary  = {'blast':'blast2'}
    if aligner_cmd in binary:
        bin_cmd = binary[aligner_cmd]
    else:
        bin_cmd = None # create_runer will know how to do

    run_align_for_seq = create_runner(kind=aligner_cmd, bin_=bin_cmd,
                           parameters=cmd_parameters)

    def _filter(sequence):
        'Giving a sequence it returns true or False depending on the exonerate'

        source_result    = run_align_for_seq(sequence)[0]
        results          = parser(source_result)
        filtered_results = FilteredAlignmentResults(match_filters=match_filters,
                                                  result_filters=result_filters,
                                                   results=results)
        try:
            #only one sequence -> only one result
            filtered_results.next()
        except StopIteration:
            #there was no result for this sequence
            return False
        return True
    return _filter

def create_length_filter(length, count_masked=True):
    '''It return a function that can check if the sequence is long enough '''

    def filter_by_length(sequence):
        'It returns true if the sequence is longer than the length'
        if count_masked:
            sequence_len = len(sequence)
        else:
            sequence_len = _count_non_masked(sequence)

        if sequence_len > length:
            return True
        else:
            return False
    return filter_by_length

def _count_non_masked(sequence):
    'It returns the unmasked seq'
    sequence_length = 0
    for letter in sequence:
        if str(letter) in  ('A', 'T', 'C', 'G'):
            sequence_length += 1
    return sequence_length

def create_adaptor_matches_filter(adaptors, number=3):
    '''It creates a filter that return Fase if the sequence have more than
    number or equal tiems.

    Adaptors can be a fasta file or a list of sequences
    '''
    # The adaptators is a file or is just a list of sequences?
    properties = dir(adaptors)
    if 'name' in properties and 'close' in properties:
        fhand = adaptors
    else:
        fhand = temp_multi_fasta_file(adaptors)
    fname = fhand.name

    parameters     =  {'target': fname}
    match_filters  = [{'kind'          : 'min_scores',
                      'score_key'      : 'similarity',
                      'min_score_value': 96},
                      {'kind'          : 'min_length',
                       'min_length_subject_%' :90 }]
    #the adaptors should be found 3 or more times
    result_filters = [{'kind': 'max_num_match_parts', 'value': 2}]

    match_filter = create_aligner_filter(aligner_cmd='exonerate',
                                    cmd_parameters=parameters,
                                    match_filters=match_filters,
                                    result_filters=result_filters )

    def filter_by_adaptor_matches(sequence):
        'It return Fase if the sequence have more than number or equal tiems'
        fhand
        return match_filter(sequence)
    return filter_by_adaptor_matches



def create_comtaminant_filter(contaminant_db):
    '''It creates a filter that return False if the sequence has a strong match
     with the database
    '''
    # This filter are bases in seqclean defaults
    parameters     =  {'database':contaminant_db, 'program':'blastn' }
    match_filters  = [{'kind'          : 'min_scores',
                      'score_key'      : 'similarity',
                      'min_score_value': 96},
                      {'kind'          : 'min_length',
                       'min_length_query_%' :60 }]
    result_filters = [{'kind': 'max_num_matches', 'value':0}]
    match_filter = create_aligner_filter(aligner_cmd='blast',
                                    cmd_parameters=parameters,
                                    match_filters=match_filters,
                                    result_filters=result_filters )

    def filter_(sequence):
        'It return Fase if the sequence have more than number or equal tiems'
        return match_filter(sequence)
    return filter_





