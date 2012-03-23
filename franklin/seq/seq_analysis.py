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
from franklin.utils.misc_utils import get_fhand
from franklin.seq.writers import temp_fasta_file
from franklin.seq.alignment_result import (filter_alignments,
                                           get_alignment_parser, BlastParser)

def get_orthologs(blast1_fhand, blast2_fhand, sub1_def_as_acc=None,
                  sub2_def_as_acc=None):
    '''It return orthologs from two pools. It needs the xml output blast of the
    pools'''
    # First we have to get hist from the first blast. We will put the in a set
    blast1_hits = set()
    for hits in get_hit_pairs_from_blast(get_fhand(blast1_fhand),
                                        sub_def_as_acc=sub1_def_as_acc):
        blast1_hits.add(hits)

    # Know we will see if the hits in the second blast in the first too
    for hits in get_hit_pairs_from_blast(get_fhand(blast2_fhand),
                                        sub_def_as_acc=sub2_def_as_acc):
        hits = (hits[1], hits[0])
        if hits in blast1_hits:
            yield hits

def get_hit_pairs_from_blast(blast_fhand, sub_def_as_acc=None, filters=None):
    '''It returns an iterator with query-subject tuples of the hits
    in the blast
    '''

    blasts = BlastParser(fhand=blast_fhand, subj_def_as_accesion=sub_def_as_acc)
    if filters is None:
        filters = [{'kind'           : 'best_scores',
                    'score_key'      : 'expect',
                    'max_score'      : 1e-20,
                    'score_tolerance': 10}]
    filtered_results = filter_alignments(blasts, config=filters)

    get_id = lambda x : x.split()[0]

    for match in filtered_results:
        try:
            query = match['query'].id
        except AttributeError:
            query = match['query'].name
        query = get_id(query)
        for match_hit in match['matches']:
            try:
                subject = match_hit['subject'].id
            except AttributeError:
                subject = match_hit['subject'].name
            subject = get_id(subject)
            yield(query, subject)

def infer_introns_for_cdna(sequence, genomic_db, genomic_seqs_index=None,
                           similar_sequence=None):
    'It infers the intron location in the cdna using est2genome'

    if not similar_sequence:
        #first we want to know where is the most similar seq in the genomic_db
        #this will speed up things
        similar_seqs = look_for_similar_sequences(sequence, database=genomic_db,
                                                  blast_program='blastn')
        if not similar_seqs:
            return []
        similar_seq = similar_seqs[0]
    else:
        similar_seq = similar_sequence
    start = similar_seq['subject_start']
    end = similar_seq['subject_end']
    try:
        similar_seq = genomic_seqs_index[similar_seq['name']]
    except KeyError:
        msg = 'Sequence %s was not found' % similar_seq['name']
        raise KeyError(msg)

    #now we run est2genome for this cdna
    cdna_file = temp_fasta_file(seqs=[sequence])
    similar_seq_file = temp_fasta_file(seqs=[similar_seq])

    #we run est2genome
    cmd = ['est2genome', cdna_file.name, similar_seq_file.name,
           '-sbegin2', str(start), '-send2', str(end), '-stdout', '-auto']
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

def look_for_similar_sequences(sequence, database, blast_program, filters=None):
    '''It return a list with the similar sequences in the database.
    First it makes a blast of the sequences against the database and then
    in joins the similar sequences
    '''
    parameters = {'database': database}

    blast_runner = create_runner(tool=blast_program, parameters=parameters)
    blast_fhand  = blast_runner(sequence)[blast_program]
    return similar_sequences_for_blast(blast_fhand, filters=filters)

def similar_sequences_for_blast(blast_fhand, filters=None):
    'It look for similar sequences in a blast result'
    #now we parse the blast
    blast_parser = get_alignment_parser('blast+')
    blast_result = blast_parser(blast_fhand)

    # We filter the results with appropiate  filters
    if filters is None:
        filters = [{'kind'     : 'score_threshold',
                    'score_key': 'similarity',
                    'min_score': 90,
                   },
                   {'kind'            : 'min_length',
                    'min_num_residues': 100,
                    'length_in_query' : True
                   }
                  ]
    alignments = filter_alignments(blast_result, config=filters)
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

def do_transitive_clustering(similar_pairs):
    'Given an iterator of similar pairs it returns the transtive clusters'

    # There are much more efficient algorithms like
    # http://en.wikipedia.org/wiki/Floyd%E2%80%93Warshall_algorithm
    clusters = []
    for pair in similar_pairs:
        # which clusters share a member with this pair?
        clusters_to_join = []
        for cluster_index in range(len(clusters)):
            cluster = clusters[cluster_index]
            if pair[0] in cluster or pair[1] in cluster:
                clusters_to_join.append(cluster_index)
        if len(clusters_to_join) == 1:
            clusters[clusters_to_join[0]].update(pair)
        elif len(clusters_to_join) == 2:
            # we join all clusters that share something
            joined_cluster = set()
            for cluster_index in clusters_to_join:
                joined_cluster.update(clusters[cluster_index])
            joined_cluster.update(pair)
            clusters_to_join.sort(reverse=True)
            for cluster_index in clusters_to_join:
                del clusters[cluster_index]
            clusters.append(joined_cluster)
        elif len(clusters_to_join) > 2:
            raise RuntimeError('More than two clusters to join. Fix me!')
        else:
            clusters.append(set(pair))
    return clusters


def do_transitive_clustering_on_blast(blast_fhand, filters=None):
    '''It does a transtive clustering given a xml blast result.

    It will look for pairs of similar in the blast. Then it will create
    clusters of sequences using the transitive property in the pairs.
    e.g. a is similar to b. b to c => a, b and c belong to the same cluster.
    '''
    similar_pairs = get_hit_pairs_from_blast(blast_fhand, filters=filters)
    return do_transitive_clustering(similar_pairs)


def select_longer_sequence_from_cluster(seqs_index, clusters):
    '''It selects one sequence for every cluster and creates an output fasta.

    The sequences should be given as a BioPython sequence index (dict like).
    '''

    # It returns the index of the maximum value of a list.
    # If there is more than one max value it will return only the first one.
    index_of_max = lambda list_: list_.index(max(list_))

    for cluster in clusters:
        cluster = tuple(cluster)
        lengths = [len(seqs_index[seq_name]) for seq_name in cluster]
        yield seqs_index[cluster[index_of_max(lengths)]]
