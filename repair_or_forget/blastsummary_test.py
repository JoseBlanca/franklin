'''
Created on 2009 mai 7

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

import unittest
from franklin.blast_summary import BlastSummaries, summarie_to_gff3
import franklin
import os.path
DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

#pylint: disable-msg=R0904
class BlastSummariesTests(unittest.TestCase):
    '''Blast Summaries tests '''
    @staticmethod
    def test_summarizer_init():
        ''' Basic init test'''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        BlastSummaries(open(fname,'r'))

    @staticmethod
    def test_no_filter():
        '''Test if no filter '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(open(fname,'r'))
        expected  = {'cCL1Contig2':3, 'cCL1Contig3':1, 'cCL1Contig4':5,
                     'cCL1Contig5':8}
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)
    @staticmethod
    def test_filter_best_expects():
        ''' Test Filter: Best expected  '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(open(fname,'r'))
        expected  = {'cCL1Contig2':2, 'cCL1Contig3':1, 'cCL1Contig4':1,
                   'cCL1Contig5':2}
        summaries.add_filter_best_expects(min_expect = 1e-4,
                                           expect_tolerance = 10)
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)

    @staticmethod
    def test_filter_expect_threshold():
        ''' Test Filter: expected treshold  '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(open(fname,'r'))
        expected  = {'cCL1Contig2':2, 'cCL1Contig4':2, 'cCL1Contig5':2}
        summaries.add_filter_expect_threshold(min_expect = 1e-34)
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)

    @staticmethod
    #pylint: disable-msg=C0103
    def test_filter_similarity_threshold():
        ''' Test Filter: similarity treshold  '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(open(fname,'r'))
        expected  = {'cCL1Contig4':1}
        summaries.add_filter_similarity_threshold(min_length = 430,
                                              min_similarity = 90)
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)
    @staticmethod
    #pylint: disable-msg=C0103
    def test_filter_compatibility_threshold():
        ''' Test Filter: filter compatibility treshold  '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(open(fname,'r'))
        expected  = {'cCL1Contig4':1}
        summaries.add_filter_compatibility_threshold(min_compatibility = 400,
                                                 max_incompatibility = 50,
                                                 min_similarity = 60)
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)
    @staticmethod
    def test_reverse_alignments():
        '''It tests that the aligment subject are in reverse strand
        And we change them'''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(open(fname,'r'))
        for summarie in summaries:
            for hit in summarie._hits:
                for hsp in hit['hsps']:
                    assert hsp['query_start'] < hsp['query_end']
                    assert hsp['subject_start'] < hsp['subject_end']

#class BlastSummaries2gff3Tests(unittest.TestCase):
#    '''It test function related to blast2gff '''
#
#    @staticmethod
#    def test_summarie_to_gff3():
#        '''It test summarie_to_gff3 '''
#        fname = os.path.join(DATA_DIR, 'blast.xml')
#        fname = '/home/peio/work_in/blast2go/SGN-U576037-Alignment.xml'
#        summaries = BlastSummaries(open(fname,'r'))
#        for summarie in summaries:
#            assert summarie_to_gff3(summarie, 'cluster', query_db="SGD",
#                                   query_regex="SGN-(\w+)",
#                                   subject_regex="gi\|(\w+)*", subject_db ="GB")

if __name__ == "__main__":
    unittest.main()
