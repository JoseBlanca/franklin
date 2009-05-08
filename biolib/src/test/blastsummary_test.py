'''
Created on 2009 mai 7

@author: peio
'''
import unittest
from biolib.blast_summary import BlastSummaries
import biolib
import os.path
DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

#pylint: disable-msg=R0904
class BlastSummariesTests(unittest.TestCase):
    '''Blast Summaries tests '''
    @staticmethod
    def test_summarizer_init():
        ''' Basic init test'''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        BlastSummaries(blast=fname)

    def test_missing_blast_file_error(self):
        ''' Missing blast file test'''
        self.failUnlessRaises(OSError, BlastSummaries, blast='foo' )

    @staticmethod
    def test_no_filter():
        '''Test if no filter '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(blast=fname)
        expected  = {'cCL1Contig2':3, 'cCL1Contig3':1, 'cCL1Contig4':5,
                     'cCL1Contig5':8}
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)
    @staticmethod
    def test_filter_best_expects():
        ''' Test Filter: Best expected  '''
        fname = os.path.join(DATA_DIR, 'blast.xml')
        summaries = BlastSummaries(blast=fname)
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
        summaries = BlastSummaries(blast=fname)
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
        summaries = BlastSummaries(blast=fname)
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
        summaries = BlastSummaries(blast=fname)
        expected  = {'cCL1Contig4':1}
        summaries.add_filter_compatibility_threshold(min_compatibility = 400,
                                                 max_incompatibility = 50,
                                                 min_similarity = 60)
        for result in summaries:
            if result.query_name in expected:
                assert expected[result.query_name] == len(result.hits)


if __name__ == "__main__":
    unittest.main()