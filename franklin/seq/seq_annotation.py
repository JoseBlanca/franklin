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

from xml.parsers.expat import ExpatError
from Bio import SeqIO
from Bio.SeqFeature import  FeatureLocation
from Bio.Alphabet import generic_dna, generic_protein

from franklin.seq.alignment_result import (BlastParser, filter_alignments,
                                           build_relations_from_aligment)
from franklin.snv.snv_annotation import (INVARIANT, SNP, DELETION, INSERTION,
                                         SNV_TYPES)
from franklin.utils.cmd_utils import  create_runner
from franklin.seq.seqs import SeqFeature, get_seq_name, Seq
from franklin.utils.seqio_utils import get_content_from_fasta
from franklin.seq.seq_analysis import infer_introns_for_cdna, get_orthologs
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.misc_utils import get_fhand
from franklin.coordsystem import CoordSystem
from tempfile import NamedTemporaryFile

def create_ortholog_annotator(blast, reverse_blast, species):
    '''It creates a function factory that calculates all the orthologs between
     crossed species. First it calculates all the orthologs'''
    blast_fhand = blast['blast']
    blast_subj_def_as_acc = blast.get('subj_def_as_acc', None)
    reverse_blast_fhand = reverse_blast['blast']
    reverse_blast_subj_def_as_acc = reverse_blast.get('subj_def_as_acc', None)

    #we index the orthologs by the first one
    orthologs = {}
    for ortholog in get_orthologs(blast_fhand, reverse_blast_fhand,
                                  blast_subj_def_as_acc,
                                  reverse_blast_subj_def_as_acc):
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
                'max_score'      : 1e-20,
                'score_tolerance': 10}]
    for blast in blasts:
        blast_fhand = blast['blast']
        if 'modifier' in blast:
            modifier = blast['modifier']
        else:
            modifier = None
        blast_fhand = get_fhand(blast_fhand)
        blast = BlastParser(fhand=blast_fhand)
        filtered_results = filter_alignments(blast, config=filters)

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

#def create_polia_annotator():
#    'It creates a function that annotates polia'
#    parameters = {'min_score':'10', 'end':'x', 'incremental_dist':'20',
#                  'fixed_dist':None}
#    runner = create_runner(tool='trimpoly', parameters=parameters)
#
#    def annotate_polia(sequence):
#        'It annotates the polia'
#        if sequence is None:
#            return
#        trimpoly_out_fhand = runner(sequence)['sequence']
#        trimp_data = trimpoly_out_fhand.readline().split()
#        end5 = int(trimp_data[2]) - 1
#        end3 = int(trimp_data[3])
#
#        if end3 != 0 or end5 != len(sequence):
#            polia = SeqFeature(location=FeatureLocation(end5, end3),
#                             type='polia', qualifiers={})
#            sequence.features.append(polia)
#        return sequence
#
#    return annotate_polia

def create_prot_change_annotator():
    '''It creates a function that caracterizes if the annotated snv produces a
    change in the protein'''
    def prot_change_annotator(sequence):
        'It annotates the change in the protein by the snv'
        if sequence is None:
            return

        for orf in sequence.get_features(kind='orf'):
            for snv in sequence.get_features(kind='snv'):
                protein_change = protein_change_annotator(sequence, orf, snv)
                if protein_change is not None:
                    snv.qualifiers['protein_change'] = protein_change
        return sequence
    return prot_change_annotator

def  protein_change_annotator(sequence, orf, snv):
    '''It annotates the protein changes and stores it in the snv qualifiers
    structure'''
    codons, location, snv_pos_in_orf = _get_codons(sequence, orf, snv)
    if codons is None:
        return None
    protein_changes = {}
    protein_changes['alleles'] = {}
    for allele in snv.qualifiers['alleles']:
        protein_changes['alleles'][allele] = {}
        allele_kind = allele[1]
        aa = codons[allele].translate()
        #print aa, codons[allele]
        protein_changes['alleles'][allele]['aa'] = aa
        if aa == '*':
            allele_type = 'stop'
        elif allele_kind in (SNP, INVARIANT):
            allele_type = 'unknown'
        elif allele_kind in  (INSERTION, DELETION):
            insert_len = len(allele[0])
            codon_pos = insert_len % 3
            if codon_pos == 0:
                allele_type = SNV_TYPES[allele_kind]
            else:
                allele_type = 'frame_change'
        protein_changes['alleles'][allele]['kind'] = allele_type

    change_kind = _get_prot_change_kind(protein_changes)
    location = _get_change_locations(snv_pos_in_orf, location, change_kind)

    protein_changes['kind'] = change_kind
    protein_changes['location'] = location
    return protein_changes

def _get_prot_change_kind(protein_change):
    'it merges the protein change kind in a general'
    kinds = [allele['kind'] for allele in protein_change['alleles'].values()]
    if 'insertion' in kinds or 'deletion' in kinds:
        kind = 'indel'
    elif 'stop' in kinds or 'frame_change' in kinds:
        kind = 'breakage'
    else:
        aas = [allele['aa'] for allele in protein_change['alleles'].values()]
        if len(set(aas)) == 1:
            kind = 'synonym'
        else:
            kind = 'substitution'
    return kind

def _get_change_locations(snv_pos_in_orf, location, kind):
    'It locates the protein change in the sequence'
    if kind  == 'indel':
        location = 'utr'
    elif location == 'orf':
        location = 'codon_%d' % ((snv_pos_in_orf + 1) % 3)
    elif location is None:
        location = 'unknown'
    return location

def _locate_codons_in_orf(sequence, orf, snv):
    'It locates the snv in the orf coordinate system'
    query_name = sequence.name
    orf_seq = orf.qualifiers['dna']
    subject_name = 'subject'
    subject_fhand = NamedTemporaryFile(suffix='.fasta')
    subject_fhand.write('>%s\n%s\n' % (subject_name, orf_seq))
    subject_fhand.flush()
    parameters   = {'subject':subject_fhand.name}
    aligner      = create_runner(tool='water', parameters=parameters,
                                 add_ext_dir=False)
    result_fhand = aligner(sequence)['water']
    relations = build_relations_from_aligment(result_fhand,
                                              query_name=sequence.name,
                                              subject_name=subject_name)
    #print relations
    coord = CoordSystem(relations=[relations])

    # snv .positions
    snv_pos = snv.location.start.position


    try:
        snv_in_orf = coord.transform(from_mol=query_name, to_mol=subject_name,
                                     position=snv_pos)
    except RuntimeError:
        snv_in_orf = None

    #print snv_in_orf_start, snv_in_orf_end

    orf_start = 0
    orf_end   = len(orf.qualifiers['dna']) - 1
    orf_start_limit_in_seq = coord.transform(from_mol=subject_name,
                                             to_mol=query_name,
                                             position=orf_start)
    orf_end_limit_in_seq   = coord.transform(from_mol=subject_name,
                                             to_mol=query_name,
                                             position=orf_end)

    if snv_in_orf is None:
        # it can be utr3, utr5, or None
        if snv_pos < orf_start_limit_in_seq:
            position = 'utr5'
        elif snv_pos > orf_end_limit_in_seq:
            position = 'in utr3'
        else:
            position = None
        codon_start = None
        snv_in_orf  = None
    else:
        start_codon_pos = snv_in_orf % 3
        codon_start     = snv_in_orf - start_codon_pos
        position        = 'orf'

    return (position, codon_start, snv_in_orf)

def _get_codons_with_alleles(codon_pos, snv_pos, alleles, sequence):
    'It get codons sequence giving section'
    codons = {}
    for allele, kind in alleles:
        #print codon_pos, snv_pos
        if kind in (SNP, INVARIANT):
            left  = sequence[codon_pos: snv_pos]
            rigth = sequence[snv_pos + 1: codon_pos + 3]
            codon = left + allele + rigth

        if kind == INSERTION:
            left  = sequence[codon_pos: snv_pos]
            rigth = sequence[snv_pos: codon_pos + 3]
            codon = left + allele + rigth

            index = snv_pos + 1
            while len(codon) % 3 != 0:
                codon += sequence[index]
                index += 1

        elif kind == DELETION :
            left  = sequence[codon_pos: snv_pos]
            del_allele = sequence[snv_pos + len(allele):snv_pos + len(allele) + len(allele)]
            codon = left + del_allele

            index = snv_pos + 1
            while len(codon) % 3 != 0:
                codon += sequence[index]
                index += 1

        codons[(allele, kind)] = codon

    return codons

def _get_codons(sequence, orf, snv):
    'It gets the codons afected by the snv'
    location, codon_pos, snv_pos = _locate_codons_in_orf(sequence, orf, snv)
    orf_seq = orf.qualifiers['dna']
    if codon_pos is not None and snv_pos is not None:
        alleles = snv.qualifiers['alleles'].keys()
        codons  = _get_codons_with_alleles(codon_pos, snv_pos, alleles, orf_seq)
    else:
        codons =  None
    return codons, location, snv_pos

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
        if start < end:
            qualifiers['strand'] = 'forward'
        else:
            qualifiers['strand'] = 'reverse'
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
    def annotate_intron(sequence):
        'It adds the orf to the SeqFeatures'
        if sequence is None:
            return
        try:
            introns = infer_introns_for_cdna(sequence=sequence,
                                          genomic_db=genomic_db,
                                          genomic_seqs_index=genomic_seqs_index)
        except KeyError as error:
            error = str(error).lstrip('u').strip("'")
            if 'not found' in error:
                error += ' in seq file %s, but present in blast db %s' % \
                                           (genomic_seqs_fhand.name, genomic_db)
            raise RuntimeError(error)

        for intron_pos in introns:
            feature = SeqFeature(location=FeatureLocation(intron_pos,
                                                          intron_pos),
                                 type='intron',
                                 qualifiers={'genomic_db':genomic_db})
            sequence.features.append(feature)
        return sequence
    return annotate_intron

def create_go_annotator(annot_fpath):
    'It annotates GOs using blast2go4pipe'

    go_annotations = _parse_b2g_output(open(annot_fpath))

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
