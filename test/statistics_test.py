'''Some tests for the sequence and contig statistics calculations.
Created on 10/07/2009
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
from StringIO import StringIO
import os
from tempfile import NamedTemporaryFile

from franklin.utils.misc_utils import float_lists_are_equal
from franklin.utils.seqio_utils import seqs_in_file
from franklin.statistics import (seq_distrib, general_seq_statistics,
                               seq_distrib_diff, histogram)
import franklin

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

def _read_distrib_file(fhand):
    '''Given an fhand with a distribution/histogram in it it returns the
    distrib and the bin_edges'''
    val_line = fhand.readline()
    distri = [float(item) for item in val_line.split(':')[1].strip().split(',')]
    bin_line = fhand.readline()
    edges = [float(item) for item in bin_line.split(':')[1].strip().split(',')]
    return distri, edges

def _check_distrib_file(fhand, distrib, bin_edges):
    'Given an fhand, a distrib list and the bin_edges, it checks they match'
    fhand.seek(0)
    fdistrib, fbin_edges = _read_distrib_file(fhand)
    assert distrib == fdistrib
    float_lists_are_equal(fbin_edges, bin_edges)

def _file_length(fhand):
    'Given an fhand it returns the file length in bytes'
    return os.stat(fhand.name)[6]

class StatisticsTest(unittest.TestCase):
    'It tests different statistics for the sequences'

    @staticmethod
    def test_sequence_length_distrib():
        'It tests that we get the correct sequence length distribution'
        fhand = StringIO('>h\nATCTCAT\n>o\nACTAGG\n>l\AGCTAGCGTAGT\n>a\nGTAT\n')
        distrib_fhand = NamedTemporaryFile(suffix='.txt')
        fig_fhand     = NamedTemporaryFile(suffix='.png')
        seqs = seqs_in_file(fhand)
        seq_distrib('seq_length_distrib', sequences=seqs,
                    distrib_fhand=distrib_fhand, plot_fhand=fig_fhand)
        #now we check the results
        distrib = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
        bin_edges = [0.0, 0.35, 0.7, 1.05, 1.4, 1.75, 2.1, 2.45, 2.8, 3.15,
                     3.5, 3.85, 4.2, 4.55, 4.9, 5.25, 5.6, 5.95, 6.3, 6.65, 7.0]
        _check_distrib_file(distrib_fhand, distrib, bin_edges)
        assert _file_length(fig_fhand) > 100

    @staticmethod
    def test_masked_seq_length_distrib():
        'It tests that we get the correct sequence length distribution'
        fhand = StringIO('>h\nATCTcat\n>o\nactagg\n>l\nActagcgtAGT\n>a\nGTAT\n')
        distrib_fhand = NamedTemporaryFile(suffix='.txt')
        fig_fhand = NamedTemporaryFile(suffix='.png')
        seqs = seqs_in_file(fhand)
        seq_distrib('masked_seq_distrib', sequences=seqs,
                    distrib_fhand=distrib_fhand, plot_fhand=fig_fhand)
        #now we check the results
        distrib = [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]
        bin_edges = [0., 0.35, 0.7, 1.05, 1.4, 1.75, 2.1, 2.45, 2.8, 3.15, 3.5,
                     3.85, 4.2, 4.55, 4.9, 5.25, 5.6, 5.95, 6.3, 6.65, 7.]
        _check_distrib_file(distrib_fhand, distrib, bin_edges)
        assert _file_length(fig_fhand) > 100

    @staticmethod
    def test_seq_quality_distrib():
        'It tests that we get the correct sequence quality distribution'
        fhand_seq = StringIO('>h\nACTG\n>o\nACTG\n>l\nACTG\n>a\nACG\n')
        fhand = StringIO('>h\n1 2 3 4 \n>o\n1 2 2 3\n>l\n1 2 3 3\n>a\n1 1 6\n')
        distrib_fhand = NamedTemporaryFile(suffix='.txt')
        fig_fhand = NamedTemporaryFile(suffix='.png')
        seqs = seqs_in_file(fhand_seq, fhand)
        seq_distrib('qual_distrib', sequences=seqs, distrib_fhand=distrib_fhand,
                    plot_fhand=fig_fhand)
        #now we check the results
        distrib = [5, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
        bin_edges = [1., 1.25, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5,
                     3.75, 4., 4.25, 4.5, 4.75, 5., 5.25, 5.5, 5.75, 6.]
        _check_distrib_file(distrib_fhand, distrib, bin_edges)
        assert _file_length(fig_fhand) > 100


    @staticmethod
    def test_general_seq_stats():
        "it test length_statistics"
        fhand_seq = StringIO('>h\nACTG\n>o\nACtG\n>l\nACtG\n>a\nAcG\n')
        fhand = StringIO('>h\n1 2 3 4 \n>o\n1 2 2 3\n>l\n1 2 3 3\n>a\n1 1 6\n')
        seqs  = seqs_in_file(fhand_seq, fhand)
        stats = general_seq_statistics(seqs)

        assert stats['seq_length']          == 15
        assert stats['seq_length_average']  == 3.75
        qualities = [1, 2, 3, 4, 1, 2, 2, 3, 1, 2, 3, 3, 1, 1, 6]
        mean_qual = sum(qualities) / float(len(qualities))
        assert stats['mean_quality']        == mean_qual
        assert stats['num_sequences']       == 4
        assert stats['max_seq_length']      == 4
        assert stats['min_seq_length']      == 3
        assert stats['seq_length_variance'] == 0.1875

    @staticmethod
    def test_seq_distrib_diff():
        'It test seq_distrib_diff'
        fhand_seq1 = StringIO('>h\nACTG\n>o\nACtG\n>l\nACtG\n>a\nAcG\n')
        fhand1 = StringIO('>h\n1 2 3 4 \n>o\n1 2 2 3\n>l\n1 2 3 3\n>a\n1 1 6\n')
        fhand_seq2 = StringIO('>h\nACTG\n>o\nACtG\n>l\nACtG\n>a\nAcG\n')
        fhand2 = StringIO('>h\n1 2 3 4 \n>o\n1 2 2 3\n>l\n1 2 3 3\n>a\n1 1 6\n')
        seqs1 = seqs_in_file(fhand_seq1, fhand1)
        seqs2 = seqs_in_file(fhand_seq2, fhand2)

        distri = seq_distrib_diff(seqs1, seqs2, 'qual_distrib')
        for num in distri['distrib']:
            assert num == 0

class HistogramTest(unittest.TestCase):
    'It checks our histogram/distribution implementation'
    @staticmethod
    def test_histogram():
        'It test our histogram implementation'
        numbers = [0, 1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10]
        bins    = 5
        numpy_distrib = ([2, 2, 4, 2, 3], [  0.,   2.,   4.,   6.,   8.,  10.])
        our_distrib = histogram(numbers, bins=bins)

        for num1, num2 in zip(numpy_distrib[0], our_distrib[0]):
            assert num1 == num2
        for num1, num2 in zip(numpy_distrib[1], our_distrib[1]):
            assert num1 == num2

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_sequence_length_distrib']
    unittest.main()
