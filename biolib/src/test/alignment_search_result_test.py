'''It tests the representation of the results from programs like blast,
ssaha2, etc. that align one sequence against a database.'''

import unittest
import os, math
import biolib
from biolib.alignment_search_result import (BlastParser,
                                            FilteredAlignmentResults,
                                            generate_score_distribution)

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

def _floats_are_equal(num1, num2):
    'Given two numbers it returns True if they are similar'
    log1 = math.log(num1)
    log2 = math.log(num2)
    return abs(log1 - log2) < 0.01

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
        assert _floats_are_equal(match_part['scores'][key],
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
                    assert _floats_are_equal(bl_match['scores'][key],
                                             expt_match['scores'][key])

class BlastParserTest(unittest.TestCase):
    'It test the blast parser'

    @staticmethod
    def test_blast_parser():
        'It test the blast parser'
        blast_file = open(os.path.join(DATA_DIR, 'blast.xml'))
        parser = BlastParser(fhand=blast_file)

        expected_results = [
            {'query':{'name':'lcl|2_0', 'description':'cCL1Contig2',
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
            {'query':{'name':'lcl|3_0', 'description':'cCL1Contig3',
                      'length':629},
            },
            {}, {}
        ]
        n_blasts = 0
        for index, blast in enumerate(parser):
            _check_blast(blast, expected_results[index])
            n_blasts += 1
        assert n_blasts == 4

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
        assert len(match_summary[query_name]) == expected[query_name]

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
        expected  = {'lcl|2_0':3, 'lcl|3_0':1, 'lcl|4_0':5,
                     'lcl|5_0':8}
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
        expected  = {'lcl|2_0':2, 'lcl|3_0':1, 'lcl|4_0':1,
                     'lcl|5_0':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
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
        expected  = {'lcl|2_0':2, 'lcl|3_0':0, 'lcl|4_0':2,
                     'lcl|5_0':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with similartiry
        filters = [{'kind'           : 'min_scores',
                    'score_key'      : 'similarity',
                    'min_score_value': 90,
                   }]
        expected  = {'lcl|2_0':0, 'lcl|3_0':0, 'lcl|4_0':1,
                     'lcl|5_0':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
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
        expected  = {'lcl|2_0':3, 'lcl|3_0':0, 'lcl|4_0':1,
                     'lcl|5_0':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with the min length given in query %
        filters = [{'kind'               : 'min_length',
                    'min_length_query_%' : 70,
                   }]
        expected  = {'lcl|2_0':0, 'lcl|3_0':0, 'lcl|4_0':2,
                     'lcl|5_0':0}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

        #with the min length given in subject %
        filters = [{'kind'                 : 'min_length',
                    'min_length_subject_%' : 0.002,
                   }]
        expected  = {'lcl|2_0':3, 'lcl|3_0':0, 'lcl|4_0':1,
                     'lcl|5_0':2}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
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
        expected  = {'lcl|2_0':0, 'lcl|3_0':0, 'lcl|4_0':1,
                     'lcl|5_0':0}
        blasts = BlastParser(fhand=blast_file)
        filtered_blasts = FilteredAlignmentResults(filters=filters,
                                                   results=blasts)
        match_summary = _summarize_matches(filtered_blasts)
        _check_match_summary(match_summary, expected)

class AlignmentSearchSimilDistribTest(unittest.TestCase):
    'It test that we can calculate the distribution of similarity'

    @staticmethod
    def test_scores_distribution():
        'We can calculate scores distributions for the alinment search result'
        #some faked test data
        result1 = {'matches':
                        [{'start':10, 'end':20,
                          'scores':{'expect':0.01},
                          'match_parts':[{'scores':{'similarity':90.0}}]
                         },
                         {'start':30, 'end':40,
                          'scores':{'expect':0.01},
                          'match_parts':[{'scores':{'similarity':80.0}}]
                         }]}
        result2 = {'matches':
                        [{'start':10, 'end':20,
                          'scores':{'expect':0.01},
                          'match_parts':[{'scores':{'similarity':60.0}}]
                         },
                        {'start':10, 'end':20,
                          'scores':{'expect':0.01},
                          'match_parts':[{'scores':{'similarity':60.1}}]
                         },
                         {'start':30, 'end':40,
                          'scores':{'expect':0.01},
                          'match_parts':[{'scores':{'similarity':80.1}}]
                         }]}
        blasts = [result1, result2]
        distrib = generate_score_distribution(results=blasts,
                                              score_key='similarity',
                                              nbins = 3)
        assert distrib['distribution'] == [22, 22, 11]
        assert distrib['bins'] == [60.0, 70.0, 80.0, 90.0]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
