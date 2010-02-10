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
from biolib.utils.cmd_utils import  create_runner
from biolib.utils.seqio_utils import get_seq_name
from biolib.seq.seqs import SeqFeature
from biolib.utils.seqio_utils import get_content_from_fasta

from Bio.SeqFeature import  FeatureLocation
from Bio.Alphabet import generic_dna, generic_protein
from Bio.Seq import Seq

def create_ortholog_annotator(blast, reverse_blast, species):
    '''It creates a function factory that calculates all the orthologs between
     crossed species. First it calculates all the orthologs'''
    blast_fhand         = blast['results']['blast']
    reverse_blast_fhand = reverse_blast['results']['blast']
    orthologs = list(_get_orthologs(blast_fhand, reverse_blast_fhand))

    def ortholog_annotator(sequence):
        'The real annotator'
        name = get_seq_name(sequence)
        seq_orthologs = []
        for ortholog in orthologs:
            if ortholog[0] == name:
                seq_orthologs.append(ortholog[1])
        sequence.annotations['%s-orthologs' % species] = seq_orthologs
        return sequence
    return ortholog_annotator

def _get_orthologs(blast1_fhand, blast2_fhand):
    '''It return orthologs from two pools. It needs the xml otput blast of the
    pools'''
    # First we have to get hist from the first blast. We will put the in a set
    blast1_hits = set()
    for hits in get_hit_pairs_fom_blast(blast1_fhand):
        blast1_hits.add(hits)

    # Know we will see if the hits in the second blast in the first too
    for hits in get_hit_pairs_fom_blast(blast2_fhand):
        hits = (hits[1], hits[0])
        if hits in blast1_hits:
            yield hits

def get_hit_pairs_fom_blast(blast1_fhand, filters=None):
    'It return a iterator with query subjetc tuples of the hist in the blast'

    blasts = BlastParser(fhand=blast1_fhand)
    if filters is None:
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

def create_description_annotator(blasts):
    '''It creates a function that return the best description from a list of
    blast results'''
    descriptions = _get_descriptions_from_blasts(blasts)
    def descrition_annotator(sequence):
        'The description annotator'
        name = get_seq_name(sequence)
        if name in descriptions:
            sequence.annotations['description'] = descriptions[name]
        return sequence
    return descrition_annotator

def _get_descriptions_from_blasts(blasts):
    '''It gets a description from a list of blast outputs.
    Blast description in the xml may be modified to remove trash. This depends
    on blast xml, so the item of the list can be a blast or a dict with the
    blast and the function to modify the description field.

    It tries to find the name in the first file, after in the second, etc'''

    seq_annot = {}
    filters = [{'kind'           : 'best_scores',
                'score_key'      : 'expect',
                'max_score_value': 1e-4,
                'score_tolerance': 10}]
    for blast in blasts:
        blast_fhand = blast['results']['blast']
        if 'modifier' in blast:
            modifier = blast['modifier']
        else:
            modifier = None

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
                if modifier is not None:
                    description = modifier(description)
                if description != "<unknown description>":
                    seq_annot[query] = description.strip()
    return seq_annot


def create_microsatelite_annotator():
    'It creates a function that'
    runner = create_runner(tool='sputnik')

    def search_ssr(sequence):
        'Do the actual search'
        srrs_out_fhand = runner(sequence)['sputnik']
        for feature in _get_features_from_sputnik(srrs_out_fhand):
            sequence.features.append(feature)
        return sequence
    return search_ssr


def _get_features_from_sputnik(fhand):
    'It parses the sputnik output'
    for line in fhand:
        line = line.strip()
        if 'score' in line:
            items = line.split()
            ssr_type = items[0]
            start = int(items[1]) - 1
            end   = int(items[3]) - 1
            score = int(items[8])
            unit  = items[10]
            yield SeqFeature(location=FeatureLocation(start, end), type='ssr',
                             qualifiers={'score':score, 'unit':unit,
                                         'type':ssr_type})
def create_orf_annotator(parameters):
    'It creates a function that annotates orfs'
    runner = create_runner(tool='estscan', parameters=parameters)

    def annotate_orf(sequence):
        'It adds the orf to the SeqFeatures'
        results = runner(sequence)
        dna_fhand = results['dna']
        prot_fhand = results['protein']
        description, seq = get_content_from_fasta(dna_fhand)[1:]
        pep = get_content_from_fasta(prot_fhand)[-1]
        start, end = description.split()[-2:]
        start, end = int(start), int(end)
        seq = Seq(seq, generic_dna)
        pep = Seq(pep, generic_protein)
        feature = SeqFeature(location=FeatureLocation(start, end), type='orf',
                             qualifiers={'dna':seq, 'pep':pep})
        sequence.features.append(feature)

    return annotate_orf






