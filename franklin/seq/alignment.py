'''
Created on 06/04/2011

@author: jose
'''

from franklin.seq.alignment_result import (get_alignment_parser,
                                           filter_alignments)
from franklin.utils.cmd_utils import create_runner
from franklin.seq.writers import temp_fasta_file

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
