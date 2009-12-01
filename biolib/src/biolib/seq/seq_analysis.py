'''
Created on 26/11/2009

@author: jose
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

from biolib.utils.cmd_utils import create_runner
from biolib.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)
from biolib.utils.collections_ import list_consecutive_pairs_iter

def _infer_introns_from_matches(alignments):
    'Given a match with several match parts it returns the introns'
    alignment = alignments.next()
    match     = alignment['matches'][0]
    hsps      = match['match_parts']

    introns = []
    direct_hsps, reverse_hsps = _separate_hsps(hsps)
    #from pprint import pprint
    #pprint (direct_hsps)
    #pprint (reverse_hsps)
    for hsps, orient in ((direct_hsps, 'direct'), (reverse_hsps, 'reverse')):
        #we sort the hsps to compare the consecutive ones
        if orient == 'direct':
            sort_by = 'query_end'
        else:
            sort_by = 'query_start'
        hsps = sorted(hsps, lambda x, y: x[sort_by] - y[sort_by])
        #pprint(hsps)
        for hsp1, hsp2 in list_consecutive_pairs_iter(hsps):
            #pprint(hsp1)
            #pprint(hsp2)
            intron = _infer_introns_form_match_parts(hsp1, hsp2)
            #print 'intron', intron
            #print '***************************'
            if intron:
                introns.append(intron)
    introns = sorted(introns, lambda x, y: x - y)
    return introns

def _infer_introns_form_match_parts(hsp1, hsp2):
    'It looks for introns between two hsps'
    blast_tolerance = 10 #in base pairs
    intron_tolerance = 0.05 #error
    #which hsp1 should start before than hsp2
    if hsp1['query_start'] > hsp2['query_start']:
        hsp1, hsp2 = hsp2, hsp1
    #we need the points that define the gap between the hsps
    # \
    #  \
    #   x 1  x 2 point 1 and 2
    #         \
    #          \
    #           \
    point1, point2 = {}, {}
    if hsp1['query_strand'] > 0:
        point1['query']   = hsp1['query_end']
        point1['subject'] = hsp1['subject_end']
        point2['query']   = hsp2['query_start']
        point2['subject'] = hsp2['subject_start']
    else:
        point1['query']   = hsp1['query_start']
        point1['subject'] = hsp1['subject_start']
        point2['query']   = hsp2['query_end']
        point2['subject'] = hsp2['subject_end']
    #print 'point1', point1
    #print 'point2', point2
    #the point 2 should be at the same position as point1 or to 3'
    # \
    #  \
    #   xooooo
    #   oooooo
    #   oooooo
    query_dif   = point2['query'] - point1['query']
    subject_dif = point2['subject'] - point1['subject']
    if query_dif + blast_tolerance < 0:
        return None
    if subject_dif + blast_tolerance < 0:
        return None
    #now the real introns
    # \
    #  \
    #   xooooo
    #     oooo
    #       oo
    #        o (this is not a straight line!)

    # this is a very strange case. when it happens there is no intron
    if float(point2['subject'] - point1['subject']) == 0:
        return None
    intron_index = (point2['subject'] - point1['subject'] - point2['query'] +
                    point1['query']) / \
                    float(point2['subject'] - point1['subject'])
    #print intron_index
    if intron_index > intron_tolerance:
        return int((point2['query'] + point1['query']) / 2.0)
    else:
        return None

def _separate_hsps(hsps):
    'It separates hsps taking into accound the query and subject strands'
    direct_hsps  = []
    reverse_hsps = []
    for hsp in hsps:
        if (hsp['query_strand'] * hsp['subject_strand']) > 0:
            direct_hsps.append(hsp)
        else:
            reverse_hsps.append(hsp)
    return direct_hsps, reverse_hsps

def infer_introns_for_cdna(sequence, genomic_db):
    '''Doing a blast with the sequences against the genomic db it infers the
    positions of introns'''

    #first we run the blast
    parameters = {'database': genomic_db, 'program':'tblastx'}

    filters = [{'kind'          : 'min_length',
                'min_length_bp' : 20}]
    blast_runner = create_runner(kind='blast', parameters=parameters)
    blast_fhand = blast_runner(sequence)[0]

    #now we parse the blast
    blast_parser = get_alignment_parser('blast')
    blast_result = blast_parser(blast_fhand)

    # We filter the results with appropiate  filters

    alignments = FilteredAlignmentResults(match_filters=filters,
                                          results=blast_result)

    #now we have to guess the introns
    introns = _infer_introns_from_matches(alignments)
    return introns

def look_for_similar_sequences(sequence, db, blast_program, filters=None):
    'It return a list with the similar sequences in the database'
    parameters = {'database': db, 'program':blast_program}

    if filters is None:
        filters = [{'kind'           : 'min_scores',
                    'score_key'      : 'similarity',
                    'min_score_value': 90,
                   },
                   {'kind'           : 'min_length',
                    'min_length_bp'  : 100,
                   }
                  ]

    blast_runner = create_runner(kind='blast', parameters=parameters)
    blast_fhand = blast_runner(sequence)[0]

    #now we parse the blast
    blast_parser = get_alignment_parser('blast')
    blast_result = blast_parser(blast_fhand)

    # We filter the results with appropiate  filters

    alignments = FilteredAlignmentResults(match_filters=filters,
                                          results=blast_result)
    try:
        alignment = alignments.next()
    except StopIteration:
        return []
    similar_seqs = []
    for match in alignment['matches']:
        similar_seqs.append(match['subject'].name)

    return similar_seqs




