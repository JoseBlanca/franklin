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
from biolib.alignment_search_result import (ExonerateParser, BlastParser,
                                            FilteredAlignmentResults)


def create_filter(aligner_cmd, cmd_parameters, match_filters=None,
                            result_filters=None):
    '''A function factory factory that creates exonerate filters.

    It returns a function that will accept a sequence and it will return
    True or False depending on the exonerate outcome.
    parameters is a dictionary and key are defined in ExonerateRunner.
    Required is only the target fasta file
    '''
    #runners = {'blast':BlastRunner, 'exonerate':ExonerateRunner}

    parsers = {'blast':BlastParser, 'exonerate':ExonerateParser}
    parser  = parsers[aligner_cmd]
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
        filtered_results = FilteredAlignmentResults(filters=match_filters,
                                                   results=results)
        try:
            result = filtered_results.next()
            if not len(result['matches']):
                return False
        except StopIteration:
            return False

        filter_result = _filtered_match_results(filters=result_filters,
                                                result=result)
        return filter_result

    return _filter


def _filtered_match_results(filters, result):
    '''It returns True or False depending if the result pass the filter or
    not'''
    if filters is None:
        return True
    num_matches = len(result['matches'])
    for filter_ in filters:
        if filter_['kind'] == 'num_matches':
            min_num_matches = filter_['value']
    if num_matches < min_num_matches:
        return False
    else:
        return True



####
