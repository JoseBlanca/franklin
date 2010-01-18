'''
This module has some functions to calculate the orthologs using just blast.
We are not going to build trees. So the result is not going to be as accurate
as it could be.

Created on 15/01/2010

@author: peio
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
from biolib.alignment_search_result import (BlastParser,
                                            FilteredAlignmentResults)


def get_orthologs(blast1_fhand, blast2_fhand):
    '''It return orthologs from two pools. It needs the xml otput blast of the
    pools'''
    # First we have to get hist from the first blast. We will put the in a set
    blast1_hits = set()
    for hits in _get_hit_pairs_fom_blast(blast1_fhand):
        blast1_hits.add(hits)

    # Know we will see if the hits in the second blast in the first too
    for hits in _get_hit_pairs_fom_blast(blast2_fhand):
        hits = (hits[1], hits[0])
        if hits in blast1_hits:
            yield hits


def _get_hit_pairs_fom_blast(blast1_fhand):
    'It return a iterator with query subjetc tuples of the hist in the blast'

    blasts = BlastParser(fhand=blast1_fhand)
    filters = [{'kind'           : 'best_scores',
                'score_key'      : 'expect',
                'max_score_value': 1e-4,
                'score_tolerance': 10}]
    filtered_results = FilteredAlignmentResults(match_filters=filters,
                                                results=blasts)
    for match in filtered_results:
        try:
            query =  match['query'].id
        except AttributeError:
            query =  match['query'].name
        for match_hit in match['matches']:
            try:
                subject =  match_hit['subject'].id
            except AttributeError:
                subject =  match_hit['subject'].name
            yield(query, subject)

def get_descriptions_from_blasts(blasts):
    '''It gets a description from a list of blast outputs.
    Blast description of the xml may be modified to remove trash. This depends
    blast xml, so the item of the list can be a blast or a dict with the blast
    and the function to modify the description field.

    It tries to find the name in the first file, after in the second, etc'''

    seq_annot = {}
    filters = [{'kind'           : 'best_scores',
                'score_key'      : 'expect',
                'max_score_value': 1e-4,
                'score_tolerance': 10}]
    for blast in blasts:
        if isinstance(blast, dict):
            blast_fhand   = blast['fhand']
            desc_modifier = blast['desc_modifier']
        else:
            blast_fhand   = blast
            desc_modifier = None
        blast = BlastParser(fhand=blast_fhand)
        filtered_results = FilteredAlignmentResults(match_filters=filters,
                                                    results=blast)
        #filtered_results = [filtered_results]
        for match in filtered_results:
            try:
                query =  match['query'].id
            except AttributeError:
                query =  match['query'].name
            if query not in seq_annot:
                match_hit = match['matches'][0]
                description = match_hit['subject'].description
                if desc_modifier is not None:
                    description = desc_modifier(description)
                if description != "<unknown description>":
                    seq_annot[query] = description
    return seq_annot












