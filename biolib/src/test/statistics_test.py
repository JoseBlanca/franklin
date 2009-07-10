'''Some tests for the sequence iters statistics calculations.
Created on 10/07/2009
'''

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

import unittest
from StringIO import StringIO
import os
from tempfile import NamedTemporaryFile

from biolib.biolib_utils import (seqs_in_file, float_lists_are_equal)
from biolib.statistics import seq_length_distrib, masked_seq_length_distrib


def read_distrib_file(fhand):
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
    fdistrib, fbin_edges = read_distrib_file(fhand)
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
        fig_fhand = NamedTemporaryFile(suffix='.png')
        seqs = seqs_in_file(fhand)
        seq_length_distrib(sequences=seqs, distrib_fhand=distrib_fhand,
                           plot_fhand=fig_fhand)
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
        fhand = StringIO('>h\nATCTcat\n>o\nactagg\n>l\AGctagcgtAGT\n>a\nGTAT\n')
        distrib_fhand = NamedTemporaryFile(suffix='.txt')
        fig_fhand = NamedTemporaryFile(suffix='.png')
        seqs = seqs_in_file(fhand)
        masked_seq_length_distrib(sequences=seqs, distrib_fhand=distrib_fhand,
                                  plot_fhand=fig_fhand)
        #now we check the results
        distrib = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1]
        bin_edges = [ 0. ,  0.3,  0.6,  0.9,  1.2,  1.5,  1.8,  2.1,  2.4,  2.7,
                      3. , 3.3,  3.6,  3.9,  4.2,  4.5,  4.8,  5.1,  5.4,  5.7,
                      6. ]
        _check_distrib_file(distrib_fhand, distrib, bin_edges)
        assert _file_length(fig_fhand) > 100


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_sequence_length_distrib']
    unittest.main()