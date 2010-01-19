'''It tests the representation of the results from programs like blast,
ssaha2, etc. that align one sequence against a database.'''

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

import biolib
from biolib.alignment_search_result import (BlastParser,
                                            FilteredAlignmentResults,
                                            alignment_results_scores,
                                            _compatible_incompatible_length,
                                             ExonerateParser)
from biolib.seq.seqs import SeqWithQuality
from biolib.utils.misc_utils import floats_are_equal

import unittest
import os
from StringIO import StringIO
from tempfile import NamedTemporaryFile
from Bio.Seq import UnknownSeq

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

def _check_sequence(sequence, expected):
    'It matches a sequence against an expected result'
    if 'name' in expected:
        assert sequence.name == expected['name']
    if 'description' in expected:
        assert sequence.description == expected['description']
    if 'length' in expected:
        assert len(sequence) == expected['length']

def _check_match_part(match_part, expected):
    'It matches a match_part against an expected result'
    assert match_part['query_start']    == expected['query_start']
    assert match_part['query_end']      == expected['query_end']
    assert match_part['query_strand']   == expected['query_strand']
    assert match_part['subject_start']  == expected['subject_start']
    assert match_part['subject_end']    == expected['subject_end']
    assert match_part['subject_strand'] == expected['subject_strand']
    for key in expected['scores']:
        assert floats_are_equal(match_part['scores'][key],
                                 expected['scores'][key])

def _check_blast(blast, expected):
    'It matches a blast results against the expected result'
    if 'query' in expected:
        _check_sequence(blast['query'], expected['query'])
    if 'matches' in expected:
        for match_index, expt_match in enumerate(expected['matches']):
            bl_match = blast['matches'][match_index]
            if 'subject' in expt_match:
                _check_sequence(bl_match['subject'],
                                expt_match['subject'])
            if 'match_parts' in expt_match:
                for match_part_index, expt_match_part in \
                                        enumerate(expt_match['match_parts']):
                    bl_match_part = bl_match['match_parts'][match_part_index]
                    _check_match_part(bl_match_part, expt_match_part)
            if 'scores' in expt_match:
                for key in expt_match['scores']:
                    assert floats_are_equal(bl_match['scores'][key],
                                             expt_match['scores'][key])

class BlastParserTest(unittest.TestCase):
    'It test the blast parser'

    @staticmethod
    def test_blast_parser():
        'It test the blast parser'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))
        parser = BlastParser(fhand=blast_file)

        expected_results = [
            {'query':{'name':'cCL1Contig2',
                      'description':"<unknown description>",
                      'length':1924},
             'matches':[
                 {'subject':{'name':'chr18',
                             'description':'No definition line found',
                             'length':19691255},
                  'scores':{'expect':4.60533e-35},
                  'match_parts':[{'query_start':276, 'query_end':484,
                                  'query_strand':-1,
                                  'subject_start':477142,
                                  'subject_end':477350,
                                  'subject_strand':1,
                                  'scores':{'expect':    4.60533e-35,
                                            'similarity':84.2,
                                            'identity':  84.2}
                                 }],
                 }
             ]
            },
            {'query':{'name':'cCL1Contig3',
                      'description':"<unknown description>",
                      'length':629},
            },
            {}, {}
        ]
        n_blasts = 0
        for index, blast in enumerate(parser):
            _check_blast(blast, expected_results[index])
            n_blasts += 1
        assert n_blasts == 4

        #with the subject id given in the xml blast
        expected_results = [
            {'query':{'name':'cCL1Contig2', 'description':'<unknown description>',
                      'length':1924}},{}, {}, {}]
        parser = BlastParser(fhand=blast_file)
        for index, blast in enumerate(parser):
            _check_blast(blast, expected_results[index])

        # Check using def as acceion in all the blasts
        #It changes depending on the blast output format. depends on version
        blast_file = open(os.path.join(DATA_DIR, 'melon_tair.xml'))
        parser = BlastParser(fhand=blast_file)
        assert parser.next()['matches'][0]['subject'].name == 'tair1'

    def test_blast_no_result(self):
        'It test that the xml output can be and empty string'
        blast_file = NamedTemporaryFile()
        blasts = BlastParser(fhand=blast_file)

        filters = [{'kind'           : 'best_scores',
                    'score_key'      : 'expect',
                    'max_score_value': 1e-4,
                    'score_tolerance': 10
                   }]
        filt_b = FilteredAlignmentResults(match_filters=filters, results=blasts)
        try:
            filt_b.next()
            self.fail
        except StopIteration:
            pass
        #match_summary = _summarize_matches(filt_b)
        #_check_match_summary(match_summary, expected)
#        match_summary = _summarize_matches(parser)
#        print match_summary


def _summarize_matches(parser):
    '''Given a alignment result parser it returns a dict with the matches for
    each query'''
    summary = {}
    for result in parser:
        query_name = result['query'].name
        matches    = result['matches']
        summary[query_name] = matches
    return summary

def _check_match_summary(match_summary, expected):
    '''Given a match summary it checks that the correct number of hits
    remain after a match filtering'''
    for query_name in expected:
        if expected[query_name] == 0:
            assert query_name not in match_summary
        else:
            assert len(match_summary[query_name]) == expected[query_name]

EXONERATE_OUTPUT = '''
" (peio@oiz: ~/work_in/tmp/exonerate)Command line: [exonerate --model affine:local query.fasta adaptators.fasta --showalignment False --showvulgar False --ryo cigar_like:%S %ql %tl\n]
Hostname: [oiz]
cigar_like:prueb2 7 31 + adaptor1 5 29 + 120 55 34 95
cigar_like:prueba 0 29 + adaptor2 5 34 + 135 54 35 95
cigar_like:prueba 27 54 + adaptor2 5 32 + 136 54 35 95
cigar_like:prueba 0 28 + adaptor3 5 32 + 126 54 35 95
cigar_like:prueba 27 54 + adaptor3 5 33 + 131 54 35 95
-- completed exonerate analysis
 '''

class ExonerateParserTest(unittest.TestCase):
    '''It tests the exonertae parser output - Cigar output '''
    @staticmethod
    def test_exonerate_parser():
        'It tests exonerate parser'
        parser = ExonerateParser(fhand=StringIO(EXONERATE_OUTPUT))

        expected_results = [
            {'query':{'name':'prueba', 'description':"<unknown description>",
                      'length':54},
             'matches':[
                 {'subject':{'name':'adaptor2','length':35},
                  'scores':{'score':136},
                  'start':0,
                  'end':29,
                  'match_parts':[{'query_start':0, 'query_end':29,
                                  'query_strand':1,
                                  'subject_start':5,
                                  'subject_end':34,
                                  'subject_strand':1,
                                  'scores':{'score':135, 'similarity':95}},
                                  {'query_start':27, 'query_end':54,
                                  'query_strand':1,
                                  'subject_start':5,
                                  'subject_end':32,
                                  'subject_strand':1,
                                  'scores':{'score':136, 'similarity':95}},
                                  ]}],
                 },
            {},{}]

        n_exonerates = 0
        for index, exonerate in enumerate(parser):
            _check_blast(exonerate, expected_results[index])
            n_exonerates += 1
        assert n_exonerates == 2




class AlignmentSearchResultFilterTest(unittest.TestCase):
    'It test that we can filter out matches from the blast or ssaha2 results'

    @staticmethod
    def test_no_filter():
        'It test the blast parser'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))
        parser = BlastParser(fhand=blast_file)
        match_summary = _summarize_matches(parser)
        #lcl|2_0 cCL1Contig2
        #lcl|3_0 cCL1Contig3
        #lcl|4_0 cCL1Contig4
        #lcl|5_0 cCL1Contig5
        expected  = {'cCL1Contig2':3, 'cCL1Contig3':1,
                     'cCL1Contig4':5, 'cCL1Contig5':8}
        _check_match_summary(match_summary, expected)

    @staticmethod
    def test_best_scores_filter():
        'We can keep the hits with the bests expects'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))
        filters = [{'kind'           : 'best_scores',
                    'score_key'      : 'expect',
                    'max_score_value': 1e-4,
                    'score_tolerance': 10
                   }]
        expected  = {'cCL1Contig2':2, 'cCL1Contig3':1,
                     'cCL1Contig4':1, 'cCL1Contig5':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

    @staticmethod
    def test_min_scores_filter():
        'We can keep the hits scores above the given one'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))

        #with evalue
        filters = [{'kind'           : 'min_scores',
                    'score_key'      : 'expect',
                    'max_score_value': 1e-34,
                   }]
        expected  = {'cCL1Contig2':2, 'cCL1Contig3':0,
                     'cCL1Contig4':2, 'cCL1Contig5':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with similartiry
        filters = [{'kind'           : 'min_scores',
                    'score_key'      : 'similarity',
                    'min_score_value': 90,
                   }]
        expected  = {'cCL1Contig2':0, 'cCL1Contig3':0,
                     'cCL1Contig4':1, 'cCL1Contig5':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

    @staticmethod
    def test_min_length_filter():
        'We can keep the hits length above the given one'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))

        #with the min length given in base pairs
        filters = [{'kind'          : 'min_length',
                    'min_length_bp' : 500,
                   }]
        expected  = {'cCL1Contig2':3, 'cCL1Contig3':0,
                     'cCL1Contig4':1, 'cCL1Contig5':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with the min length given in query %
        filters = [{'kind'               : 'min_length',
                    'min_length_query_%' : 70,
                   }]
        expected  = {'cCL1Contig2':0, 'cCL1Contig3':0,
                     'cCL1Contig4':2, 'cCL1Contig5':0}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with the min length given in subject %
        filters = [{'kind'                 : 'min_length',
                    'min_length_subject_%' : 0.002,
                   }]
        expected  = {'cCL1Contig2':3, 'cCL1Contig3':0,
                     'cCL1Contig4':1, 'cCL1Contig5':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

    @staticmethod
    def test_compatib_threshold_filter():
        'We can keep the hits compatible enough'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))
        #with the min length given in subject %
        filters = [{'kind'                : 'compatibility',
                    'min_compatibility'   : 400,
                    'max_incompatibility' : 50,
                    'min_similarity'      : 60
                   }]
        expected  = {'cCL1Contig2':0, 'cCL1Contig3':0,
                     'cCL1Contig4':1, 'cCL1Contig5':0}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

class ResultFilterTests(unittest.TestCase):
    'It tests results filters . It tests that we can filter some results'
    @staticmethod
    def test_min_matches_number_filter():
        'It tests that we can filter some results with too many matches'
        result1 = {'query':'query',
                   'matches': [{},{}]}
        result2 = {'query':None,
                   'matches': [{},{},{}]}
        results = iter([result1, result2])
        filter_ = {'kind': 'max_num_matches',
                   'value': 2}
        filtered_results = FilteredAlignmentResults(result_filters=[filter_],
                                                    results=results)
        assert list(filtered_results)[0] is result1

    @staticmethod
    def test_at_least_one_match():
        'A result with no matches left should be removed'
        match = {'subject':None,
                 'start':10, 'end':20,
                 'scores':{'expect':0.01},
                 'match_parts':[]}
        result1 = {'query':'query',
                   'matches': [match]}
        filters = [{'kind'           : 'min_scores',
                    'score_key'      : 'expect',
                    'max_score_value': 1e-34,
                   }]
        blasts = iter([result1])
        filtered_blasts = FilteredAlignmentResults(match_filters=filters,
                                                   results=blasts)
        assert len(list(filtered_blasts)) == 0


class IncompatibleLengthTest(unittest.TestCase):
    'It checks that we calculate the incompatible length ok'
    @staticmethod
    def test_incompatible_length():
        'It checks that we calculate the incompatible length ok'
        #        012345678901234567890
        #query   ---------------------
        #query             <--> <---->
        #simil                 90%
        #subject           <--> <---->
        #subject           ----------------------------------
        query   = SeqWithQuality(seq=UnknownSeq(length=21))
        subject = SeqWithQuality(seq=UnknownSeq(length=32))
        match_part1 = {'scores':{'similarity':90.0},
                       'query_start'   : 10,
                       'query_end'     : 13,
                       'query_strand'  : 1,
                       'subject_start' : 0,
                       'subject_end'   : 3,
                       'subject_strand': 1
                      }
        match_part2 = {'scores':{'similarity':90.0},
                       'query_start'   : 15,
                       'query_end'     : 20,
                       'query_strand'  : 1,
                       'subject_start' : 5,
                       'subject_end'   : 10,
                       'subject_strand': 1
                      }
        match = {'subject':subject,
                 'start':10, 'end':20,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part1, match_part2]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (10, 1)

        match = {'subject':subject,
                 'start':15, 'end':20,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part2]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (6, 5)

        match = {'subject':subject,
                 'start':10, 'end':13,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part1]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (4, 7)
        #                 0123456789012345678901
        #        0123456789012345678901
        #query            ---------------------
        #query             <-->  <--->
        #simil                 90%
        #subject           <-->  <--->
        #subject ----------------------
        query   = SeqWithQuality(seq=UnknownSeq(length=21))
        subject = SeqWithQuality(seq=UnknownSeq(length=22))
        match_part1 = {'scores':{'similarity':90.0},
                       'query_start'   : 1,
                       'query_end'     : 4,
                       'query_strand'  : 1,
                       'subject_start' : 10,
                       'subject_end'   : 13,
                       'subject_strand': 1
                      }
        match_part2 = {'scores':{'similarity':90.0},
                       'query_start'   : 7,
                       'query_end'     : 11,
                       'query_strand'  : 1,
                       'subject_start' : 16,
                       'subject_end'   : 20,
                       'subject_strand': 1
                      }
        match = {'subject':subject,
                 'start':1, 'end':11,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part1, match_part2]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (9, 4)

        match = {'subject':subject,
                 'start':1, 'end':4,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part1]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (4, 9)

        match = {'subject':subject,
                 'start':7, 'end':11,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part2]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (5, 8)

        #                 0123456789012345678901
        #        0123456789012345678901
        #query            ---------------------
        #query             <--<-->--->
        #simil                 90%
        #subject           <--<-->--->
        #subject ---------------------
        query   = SeqWithQuality(seq=UnknownSeq(length=21))
        subject = SeqWithQuality(seq=UnknownSeq(length=21))
        match_part1 = {'scores':{'similarity':90.0},
                       'query_start'   : 1,
                       'query_end'     : 7,
                       'query_strand'  : 1,
                       'subject_start' : 10,
                       'subject_end'   : 16,
                       'subject_strand': 1
                      }
        match_part2 = {'scores':{'similarity':90.0},
                       'query_start'   : 4,
                       'query_end'     : 11,
                       'query_strand'  : 1,
                       'subject_start' : 13,
                       'subject_end'   : 20,
                       'subject_strand': 1
                      }
        match = {'subject':subject,
                 'start':1, 'end':11,
                 'scores':{'expect':0.01},
                 'match_parts':[match_part1, match_part2]}

        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (11, 1)

        #query
        #      1----34  31----94
        #s   1              x
        #u   |             x
        #b   |            x
        #j  66           x
        #e
        #c 191    x
        #t   |   x
        #    |  x
        #  224 x
        query   = SeqWithQuality(seq=UnknownSeq(length=110))
        subject = SeqWithQuality(seq=UnknownSeq(length=250))
        match = {'match_parts': [{'query_strand': -1, 'subject_end': 66,
                                'subject_start': 1, 'query_start': 31,
                                'query_end': 94,
                                'scores': {'expect': 2.47e-21,
                                           'identity': 96.96,
                                           'similarity': 96.96},
                                'subject_strand': 1},
                               {'query_strand': -1, 'subject_end': 224,
                                'subject_start': 191, 'query_start': 1,
                                'query_end': 34,
                                'scores': {'expect': 8.7e-12,
                                           'identity': 100.0,
                                           'similarity': 100.0},
                                           'subject_strand': 1}],
               'start': 1, 'scores': {'expect': 2.4e-21},
               'end': 94,
               'subject': subject}
        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (64, 32)

        query   = SeqWithQuality(seq=UnknownSeq(length=180))
        subject = SeqWithQuality(seq=UnknownSeq(length=500))
        match = {'match_parts': [{'query_strand': -1, 'subject_end': 289,
                                 'subject_start': 158, 'query_start': 34,
                                 'query_end': 166,
                                 'scores': {'expect': 4.24e-65,
                                            'identity': 98.49,
                                            'similarity': 98.4},
                                 'subject_strand': 1},
                                 {'query_strand': -1, 'subject_end': 477,
                                  'subject_start': 447, 'query_start': 1,
                                  'query_end': 31,
                                  'scores': {'expect': 5.6e-10,
                                             'identity': 100.0,
                                             'similarity': 100.0},
                                             'subject_strand': 1}],
                 'start': 1, 'scores': {'expect': 4.24e-65},
                 'end': 166, 'subject': subject}
        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (133, 47)

        query   = SeqWithQuality(seq=UnknownSeq(length=260))
        subject = SeqWithQuality(seq=UnknownSeq(length=110))
        match = {'match_parts': [{'query_strand': -1, 'subject_end': 37,
                                  'subject_start': 1, 'query_start': 215,
                                  'query_end': 251,
                                  'scores': {'expect': 1.31e-13,
                                             'identity': 100.0,
                                             'similarity': 100.0},
                                  'subject_strand': 1},
                                 {'query_strand': -1, 'subject_end': 100,
                                  'subject_start': 72, 'query_start': 223,
                                  'query_end': 251,
                                  'scores': {'expect': 8.0e-09,
                                             'identity': 100.0,
                                             'similarity': 100.0},
                                  'subject_strand': 1}],
                 'start': 215, 'scores': {'expect': 1.34e-13},
                 'end': 251, 'subject': subject}
        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (37, 75)

        query   = SeqWithQuality(seq=UnknownSeq(length=200))
        subject = SeqWithQuality(seq=UnknownSeq(length=200))
        match = {'match_parts': [{'query_strand': -1, 'subject_end': 187,
                                  'subject_start': 127, 'query_start': 1,
                                  'query_end': 65,
                                  'scores': {'expect': 3e-20,
                                             'identity': 93.8,
                                             'similarity': 93.8},
                                  'subject_strand': 1},
                                 {'query_strand': -1, 'subject_end': 65,
                                  'subject_start': 5, 'query_start': 127,
                                  'query_end': 185,
                                  'scores': {'expect': 1.8e-18,
                                             'identity': 96.7,
                                             'similarity': 96.7},
                                             'subject_strand': 1}],
                 'start': 1, 'scores': {'expect': 3e-20},
                 'end': 185, 'subject': subject}
        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (124, 67)

        query   = SeqWithQuality(seq=UnknownSeq(length=220))
        subject = SeqWithQuality(seq=UnknownSeq(length=220))
        match = {'match_parts': [{'query_strand': -1, 'subject_end': 187,
                                  'subject_start': 100, 'query_start': 69,
                                  'query_end': 156,
                                  'scores': {'expect': 8.3e-44,
                                             'identity': 100.0,
                                             'similarity': 100.0},
                                  'subject_strand': 1},
                                 {'query_strand': -1, 'subject_end': 77,
                                  'subject_start': 47, 'query_start': 179,
                                  'query_end': 209,
                                  'scores': {'expect': 8.5e-10,
                                             'identity': 100.0,
                                             'similarity': 100.0},
                                             'subject_strand': 1}],
                 'start': 69, 'scores': {'expect': 8.3e-44},
                 'end': 209, 'subject': subject}
        compat, incompat = _compatible_incompatible_length(match, query)
        assert (compat, incompat) == (119, 66)

class AlignmentSearchSimilDistribTest(unittest.TestCase):
    'It test that we can calculate the distribution of similarity'

    @staticmethod
    def test_alignment_results_scores():
        'It tests that we can get all the scores from a result'
        query    = SeqWithQuality(seq=UnknownSeq(length=25), name='query')
        subject1 = SeqWithQuality(seq=UnknownSeq(length=32), name='hola')
        subject2 = SeqWithQuality(seq=UnknownSeq(length=32), name='query')
        match_part1 = {'scores':{'similarity':90.0},
                       'query_start'   : 10,
                       'query_end'     : 20,
                       'query_strand'  : 1,
                       'subject_start' : 0,
                       'subject_end'   : 10,
                       'subject_strand': 1
                      }
        match_part2 = {'scores':{'similarity':80.0},
                       'query_start'   : 10,
                       'query_end'     : 20,
                       'query_strand'  : 1,
                       'subject_start' : 0,
                       'subject_end'   : 10,
                       'subject_strand': 1
                      }
        result1 = {'query':query,
                   'matches':
                        [{'subject':subject1,
                          'start':10, 'end':20,
                          'scores':{'expect':0.01},
                          'match_parts':[match_part1]
                         },
                         {'subject':subject2,
                          'start':10, 'end':20,
                          'scores':{'expect':0.01},
                          'match_parts':[match_part2]
                         }]}
        results = [result1]
        #filtering the query-subj with the same names
        scores = alignment_results_scores(results, ['similarity'])
        assert scores == [90.0]
        #without the filtering
        scores = alignment_results_scores(results, ['similarity'],
                                               filter_same_query_subject=False)
        assert scores == [90.0, 80.0]

        #we can ask for several scores at once
        scores = alignment_results_scores(results, ['similarity', 'similarity'],
                                                filter_same_query_subject=False)
        assert scores == [[90.0, 80.0], [90.0, 80.0]]

        #we can also ask for a derived score
        scores = alignment_results_scores(results, ['similarity',
                                                    'd_incompatibility'],
                                                filter_same_query_subject=False)
        assert scores == [[90.0, 80.0], [16.0, 16.0]]



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
