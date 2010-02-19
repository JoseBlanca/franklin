'''
Created on 26/11/2009

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

from franklin.utils.cmd_utils import create_runner, call
from franklin.seq.writers import temp_fasta_file
from franklin.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)
from franklin.utils.collections_ import list_consecutive_pairs_iter
from franklin.seq.seq_annotation import get_hit_pairs_fom_blast
from tempfile import NamedTemporaryFile
def _infer_introns_from_matches(alignments):
    'Given a match with several match parts it returns the introns'
    try:
        alignment = alignments.next()
    except StopIteration:
        return []
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

def _infer_introns_for_cdna_blast(sequence, genomic_db):
    '''Doing a blast with the sequences against the genomic db it infers the
    positions of introns'''
    #first we run the blast
    parameters = {'database': genomic_db, 'program':'tblastx'}

    filters = [{'kind'          : 'min_length',
                'min_length_bp' : 20}]
    blast_runner = create_runner(tool='blast', parameters=parameters)
    blast_fhand = blast_runner(sequence)['blast']

    #now we parse the blast
    blast_parser = get_alignment_parser('blast')
    blast_result = blast_parser(blast_fhand)
    # We filter the results with appropiate  filters

    alignments = FilteredAlignmentResults(match_filters=filters,
                                          results=blast_result)
    #now we have to guess the introns
    introns = _infer_introns_from_matches(alignments)
    return introns

def _infer_introns_for_cdna_est2genome(sequence, genomic_db,
                                                            genomic_seqs_index):
    'It infers the intron location in the cdna using est2genome'

    #first we want to know where is the most similar seq in the genomic_db
    #this will speed up things
    similar_seqs = look_for_similar_sequences(sequence, database=genomic_db,
                                              blast_program='blastn')
    if not similar_seqs:
        return []
    similar_seq = similar_seqs[0]
    start = similar_seq['subject_start']
    end   = similar_seq['subject_end']
    similar_seq = genomic_seqs_index[similar_seq['name']]

    #now we run est2genome for this cdna
    cdna_file = temp_fasta_file(seqs=[sequence])
    similar_seq_file = temp_fasta_file(seqs=[similar_seq])

    #we run est2genome
    cmd = ['est2genome', cdna_file.name, similar_seq_file.name,
           '-sbegin2', str(start), '-send2', str(end), '-stdout', '-auto']
#    # Sometimes est2genome fails randomly, so we repeat the call once
#    try:
#        stdout, stderr, retcode = call(cmd)
#    except OSError:
    stdout, stderr, retcode = call(cmd)

    if retcode:
        msg = 'There was an error running est2genome: ' + stderr
        raise RuntimeError(msg)

    #parse est2genome
    result = est2genome_parser(stdout)

    #get_introns_from parser_result
    return result['cdna']['introns']

def est2genome_parser(output):
    '''It parses est2genome output'''
    result = {'cdna':{'introns':[], 'exons':[]},
              'genomic':{'introns':[], 'exons':[]}}
    for line in output.split('\n'):
        line = line.strip()
        if not line:
            continue
        if line.startswith('Exon'):
            items = line.split()
            genomic = {'start': int(items[3]), 'end': int(items[4])}
            cdna = {'start': int(items[6]), 'end': int(items[7])}
            result['cdna']['exons'].append(cdna)
            result['genomic']['exons'].append(genomic)
        elif line.startswith('-Intron') or line.startswith('+Intron'):
            items = line.split()
            genomic = {'start': int(items[3]), 'end': int(items[4])}
            cdna = result['cdna']['exons'][-1]['end']
            result['cdna']['introns'].append(cdna)
            result['genomic']['introns'].append(genomic)
    return result

def infer_introns_for_cdna(sequence, genomic_db, genomic_seqs_index=None,
                           method='est2genome'):
    '''Doing a blast with the sequences against the genomic db it infers the
    positions of introns'''
    if method == 'blast':
        print 'The blast method for looking for introns is not well tested'
        return _infer_introns_for_cdna_blast(sequence, genomic_db)
    else:
        return _infer_introns_for_cdna_est2genome(sequence, genomic_db,
                                                  genomic_seqs_index)



def look_for_similar_sequences(sequence, database, blast_program, filters=None):
    'It return a list with the similar sequences in the database'
    parameters = {'database': database, 'program':blast_program}

    blast_runner = create_runner(tool='blast', parameters=parameters)
    blast_fhand = blast_runner(sequence)['blast']
    return similar_sequences_for_blast(blast_fhand, filters=filters)

def similar_sequences_for_blast(blast_fhand, filters=None):
    'It look fro similar sequences ina blast result'
    #now we parse the blast
    blast_parser = get_alignment_parser('blast')
    blast_result = blast_parser(blast_fhand)

    # We filter the results with appropiate  filters
    if filters is None:
        filters = [{'kind'           : 'min_scores',
                    'score_key'      : 'similarity',
                    'min_score_value': 90,
                   },
                   {'kind'           : 'min_length',
                    'min_length_bp'  : 100,
                   }
                  ]
    alignments = FilteredAlignmentResults(match_filters=filters,
                                          results=blast_result)
    try:
        alignment = alignments.next()
    except StopIteration:
        return []
    similar_seqs = []
    for match in alignment['matches']:
        #to which sequence our query is similar?
        name = match['subject'].name
        similar_seqs.append({'name':name,
                             'subject_start': match['subject_start'],
                             'subject_end':   match['subject_end'],
                             'query_start':   match['start'],
                             'query_end':     match['end']
                             })
    return similar_seqs

def build_sequence_clusters(aligner_config, filters=None):
    '''It builds sequence clusters. It need sequences or the blast results of
    all against all alignment'''

    alignment_result = aligner_config['results']['blast']
    if filters is None:
        filters = [{'kind' : 'min_scores',
                    'score_key'      : 'similarity',
                    'min_score_value': 96},
                   {'kind'          : 'min_length',
                    'min_length_bp' : 15}]
    pairs = get_hit_pairs_fom_blast(alignment_result, filters=filters)
    # run of to run tclust
    input_fhand  = NamedTemporaryFile()
    for pair1, pair2 in pairs:
        if pair1 != pair2:
            input_fhand.write("%s\t%s\n" % (pair1, pair2))
    input_fhand.flush()
    #print open(input_fhand.name).read()
    cmd = ['tclust', input_fhand.name]
    tclust_result = call(cmd, raise_on_error=True)[0]
    if tclust_result:
        return _parse_tclust_result(tclust_result)
    else:
        return None

def _parse_tclust_result(result):
    '''It parses tcust result. It return a list of clusters as listt of elements
    [['element1', 'element2'], ['element3', 'element4']]
    '''
    clusters = []
    for line in result.split('\n'):
        line = line.strip()
        if not line or line[0] == '>':
            continue
        clusters.append(line.split())
    return clusters

