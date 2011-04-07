'''
Created on 06/04/2011

@author: jose
'''
from __future__ import division
from franklin.seq.alignment_result import (get_alignment_parser,
                                           filter_alignments,
                                           _fix_matches)
from franklin.utils.cmd_utils import create_runner
from franklin.seq.writers import temp_fasta_file
from franklin.seq.readers import seqs_in_file
from Bio.pairwise2 import align

def _seq_to_fasta_fhand(seq):
    'Given a fhand or Seq object it returns a fhand'
    if 'file' in seq.__class__.__name__.lower():
        return seq
    return temp_fasta_file(seqs=[seq])

class ExonerateAligner(object):
    'It aligns sequences using exonerate'
    def __init__(self, subject, parameters=None, filters=None):
        'It inits the class'

        if parameters is None:
            parameters = {}
        self._filters = filters

        self._parser  = get_alignment_parser('exonerate')

        self._subject_fhand = _seq_to_fasta_fhand(subject)
        parameters['target'] = self._subject_fhand.name
        self._aligner = create_runner(tool='exonerate', parameters=parameters)

    def do_alignment(self, query):
        'It returns an alignment with this query'

        alignment_fhand = self._aligner(query)['exonerate']
        # We need to parse the result
        alignments = self._parser(alignment_fhand)

        # We filter the results with appropriate filters
        if self._filters is not None:
            alignments = filter_alignments(alignments, config=self._filters)
        return alignments

class BlastAligner(object):
    'An aligner capable of aligning sequences using blast'
    def __init__(self, subject=None, database=None, program='blastn',
                 parameters=None, filters=None):
        '''It inits the class.

        Query should be a sequence and subject can be one or several.
        subject could be an fhand (fasta) or an string
        '''
        if subject is None and database is None:
            raise ValueError('Either subject or database should be given')
        if subject is not None and database is not None:
            msg = 'subject and database can not be given at the same time'
            raise ValueError(msg)

        if parameters is None:
            parameters = {}
        self._filters = filters

        if subject is not None:
            parameters['alig_format'] = 0
            self._parser  = get_alignment_parser('blast_text')
            self._subject_fhand = _seq_to_fasta_fhand(subject)
            parameters['subject'] = self._subject_fhand.name
        if database is not None:
            parameters['database'] = database
            parameters['alig_format'] = 5
            self._parser  = get_alignment_parser('blast')
        self._program = program
        self._aligner = create_runner(tool=program, parameters=parameters)

    def do_alignment(self, query):
        'It returns an alignment with this query'

        alignment_fhand = self._aligner(query)[self._program]
        # We need to parse the result
        alignments = self._parser(alignment_fhand)

        # We filter the results with appropriate filters
        if self._filters is not None:
            alignments = filter_alignments(alignments, config=self._filters)
        return alignments

def _seq_to_seqwithqualities(seq):
    'Given a file Seq or [Seq] return a list of strs'
    if 'file' in seq.__class__.__name__.lower():
        return list(seqs_in_file(seq))
    if not isinstance(seq, list) and not isinstance(seq, tuple):
        return [seq]

def match_parts_from_biopython_alignment(alignments, query_strand,
                                         subject_strand):
    'It returns a match_part list from a biopython alignment'
    GAP = '-'
    match_parts = []
    for alignment in alignments:
        query, subject, score, alig_start, alig_end = alignment
        alig_end -= 1 #coded in that way in biopython
        query_start   = alig_start
        subject_start = alig_start
        mismatches    = 0
        for pos in range(len(query)):
            query_char = query[pos].upper()
            subject_char = subject[pos].upper()
            if pos >= alig_end:
                pass
            elif alig_start <= pos:
                if (subject_char != query_char and
                    subject_char != GAP
                    and query_char != GAP):
                    mismatches += 1
            else:
                if query_char == GAP:
                    query_start   -= 1
                if subject_char == GAP:
                    subject_start -= 1
        query_end = query_start + alig_end
        subject_end = subject_start + alig_end
        assert query_start >= 0
        assert subject_start >= 0
        alig_len = alig_end - alig_start
        identity = (alig_len - mismatches)/alig_len * 100
        match_part = {'query_start':   query_start,
                      'query_end':     query_end,
                      'query_strand':  query_strand,
                      'subject_start': subject_start,
                      'subject_end':   subject_end,
                      'subject_strand':subject_strand,
                      'scores':{'identity': identity,
                                'score':    score}}
        match_parts.append(match_part)
    match_parts.sort(key=lambda x: x['scores']['score'], reverse=True)
    return match_parts


def sw_align(query, subject):
    'It aligns two sequences'
    query_string = str(query.seq)
    subject_string = str(subject.seq)
    reward = 2
    penalty = -1
    gapopen = -1
    gapextend = -.1
    penalize_end_gaps = False
    local = True
    function = align.localms if local else align.globalms

    alignment = function(query_string, subject_string, reward, penalty,
                        gapopen, gapextend, penalize_end_gaps=penalize_end_gaps)
    match_parts = match_parts_from_biopython_alignment(alignment,
                                                       query_strand=1,
                                                       subject_strand=1)
    alignment = {'query':query,
                 'matches':[{'subject':subject,
                             'match_parts':match_parts}]}
    _fix_matches(alignment, score_keys=['score'])
    return alignment

class SWAligner(object):
    'An aligner capable of aligning sequences using Smith Waterman'
    def __init__(self, subject=None, parameters=None, filters=None):
        '''It inits the class.

        Query should be a sequence and subject can be one or several.
        subject could be an fhand (fasta) or an string
        '''
        if parameters is None:
            parameters = {}
        self._filters = filters

        self._subjects = _seq_to_seqwithqualities(subject)

    def do_alignment(self, query):
        'It returns an alignment with this query'
        alignments = []
        for subject in self._subjects:
            alignment = sw_align(query, subject)
            alignments.append(alignment)
        if self._filters is not None:
            alignments = filter_alignments(alignments, config=self._filters)
        return alignments
