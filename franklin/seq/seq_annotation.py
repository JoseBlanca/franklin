'''
This module has some functions to calculate the orthologs using just blast.
We are not going to build trees. So the result is not going to be as accurate
as it could be.

Created on 15/01/2010

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

import os, tempfile
from xml.parsers.expat import ExpatError
from Bio import SeqIO
from Bio.SeqFeature import  FeatureLocation
from Bio.Alphabet import generic_dna, generic_protein

from franklin.alignment_search_result import (BlastParser,
                                            FilteredAlignmentResults)
from franklin.utils.cmd_utils import  create_runner, b2gpipe_runner
from franklin.seq.seqs import SeqFeature, get_seq_name, Seq
from franklin.utils.seqio_utils import get_content_from_fasta
from franklin.seq.seq_analysis import infer_introns_for_cdna, get_orthologs
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.misc_utils import get_fhand

def create_ortholog_annotator(blast, reverse_blast, species):
    '''It creates a function factory that calculates all the orthologs between
     crossed species. First it calculates all the orthologs'''
    blast_fhand = blast['blast']
    reverse_blast_fhand = reverse_blast['blast']

    #we index the orthologs by the first one
    orthologs = {}
    for ortholog in get_orthologs(blast_fhand, reverse_blast_fhand):
        if ortholog[0] not in orthologs:
            orthologs[ortholog[0]] = [ortholog[1]]
        else:
            orthologs[ortholog[0]].append(ortholog[1])

    def ortholog_annotator(sequence):
        'The real annotator'
        if sequence is None:
            return
        name = get_seq_name(sequence)
        try:
            sequence.annotations['%s-orthologs' % species] = orthologs[name]
        except KeyError:
            pass
        return sequence
    return ortholog_annotator

def create_description_annotator(blasts):
    '''It creates a function that return the best description from a list of
    blast results'''
    descriptions = _get_descriptions_from_blasts(blasts)
    def descrition_annotator(sequence):
        'The description annotator'
        if sequence is None:
            return
        name = get_seq_name(sequence)
        if name in descriptions:
            sequence.description = 'Similar to %s' % descriptions[name]
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
        blast_fhand = blast['blast']
        if 'modifier' in blast:
            modifier = blast['modifier']
        else:
            modifier = None
        blast_fhand = get_fhand(blast_fhand)
        blast = BlastParser(fhand=blast_fhand)
        filtered_results = FilteredAlignmentResults(match_filters=filters,
                                                    results=blast)
        #filtered_results = [filtered_results]
        try:
            for match in filtered_results:
                try:
                    query = match['query'].id
                except AttributeError:
                    query = match['query'].name
                if query not in seq_annot:
                    match_hit = match['matches'][0]
                    description = match_hit['subject'].description
                    if modifier is not None:
                        description = modifier(description)
                    if description != "<unknown description>":
                        seq_annot[query] = description.strip()
        except ExpatError as error:
            msg = str(error) + ':%s' % blast_fhand.name
            raise ExpatError(msg)
    return seq_annot

def create_microsatellite_annotator():
    'It creates a function that'
    runner = create_runner(tool='sputnik')

    def search_ssr(sequence):
        'Do the actual search'
        if sequence is None:
            return
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
            end = int(items[3]) - 1
            score = int(items[8])
            unit = items[10]
            yield SeqFeature(location=FeatureLocation(start, end),
                             type='microsatellite',
                             qualifiers={'score':score, 'unit':unit,
                                         'type':ssr_type})
def create_orf_annotator(parameters):
    'It creates a function that annotates orfs'
    runner = create_runner(tool='estscan', parameters=parameters)

    def annotate_orf(sequence):
        'It adds the orf to the SeqFeatures'
        if sequence is None:
            return
        results = runner(sequence)
        dna_fhand = results['dna']
        prot_fhand = results['protein']
        description, seq = get_content_from_fasta(dna_fhand)[1:]
        pep = get_content_from_fasta(prot_fhand)[-1]
        prot_fhand.close()
        dna_fhand.close()
        # If there is no description, ther is no org
        if description is None:
            return sequence

        items = description.split()
        start = int(items[1])
        end = int(items[2])
        seq = Seq(seq, generic_dna)
        pep = Seq(pep, generic_protein)
        qualifiers = {'dna':seq, 'pep':pep}
        if len(items) > 3:
            qualifiers['strand'] = 'reverse'
        else:
            qualifiers['strand'] = 'forward'
        feature = SeqFeature(location=FeatureLocation(start, end), type='orf',
                             qualifiers=qualifiers)

        sequence.features.append(feature)
        return sequence
    return annotate_orf

def create_cdna_intron_annotator(genomic_db, genomic_seqs_fhand):
    'It creates a function that annotates introns in cdna matching with genomic'
    genomic_seqs_fhand = get_fhand(genomic_seqs_fhand)
    genomic_seqs_index = SeqIO.index(genomic_seqs_fhand.name,
                                     guess_seq_file_format(genomic_seqs_fhand))
    def annotate_orf(sequence):
        'It adds the orf to the SeqFeatures'
        if sequence is None:
            return
        introns = infer_introns_for_cdna(sequence=sequence,
                                          genomic_db=genomic_db,
                                          genomic_seqs_index=genomic_seqs_index)

        for intron_pos in introns:
            feature = SeqFeature(location=FeatureLocation(intron_pos,
                                                          intron_pos),
                                 type='intron',
                                 qualifiers={'genomic_db':genomic_db})
            sequence.features.append(feature)
        return sequence
    return annotate_orf

def create_go_annotator(blast, annot_fpath=None, dat_fpath=None,
                        java_memory=None, prop_fpath=None):
    'It annotates GOs using blast2go4pipe'

    if annot_fpath is None:
        fhand, annot_fpath = tempfile.mkstemp()
        os.close(fhand)
    blast = get_fhand(blast)
    b2gpipe_runner(blast, annot_fpath=annot_fpath, dat_fpath=dat_fpath,
                   java_memory=java_memory, prop_fpath=prop_fpath)
    go_annotations = _parse_b2g_output(open(annot_fpath))
    if annot_fpath is None:
        os.remove(annot_fpath)

    def go_annotator(sequence):
        'The annotator'
        if sequence is None:
            return
        seq_name = get_seq_name(sequence)
        if  seq_name in go_annotations:
            sequence.annotations['GOs'] = go_annotations[seq_name]
        return sequence
    return go_annotator

def _parse_b2g_output(annot_fhand):
    'It parses blas2go annot file'
    annotations = {}
    for line in annot_fhand:
        line = line.strip()
        if not line:
            continue
        items = line.split('\t')
        name = items[0]
        go = items[1]
        if name not in annotations:
            annotations[name] = []
        annotations[name].append(go)
    return annotations








