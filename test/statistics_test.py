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
from franklin.statistics import CachedArray, histogram
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

class CachedArrayTest(unittest.TestCase):
    'It tests the store for iterables'
    @staticmethod
    def test_store():
        'It test the store for iterables'
        item_list = [1, 3, 4]
        storage = CachedArray(typecode='I')
        storage.extend(item_list)
        assert list(storage) == item_list
        assert list(storage) == item_list

    @staticmethod
    def test_store_to_disk():
        'It test the store for iterables saving to disk'
        item_list = [1, 2, 3]
        item_list2 = [4, 5, 6]
        storage = CachedArray(typecode='I')
        storage.extend(item_list)
        storage.to_disk()
        storage.extend(item_list2)
        assert list(storage) == [1, 2, 3, 4, 5, 6]

    @staticmethod
    def test_basic_statistics():
        'It test the max, min avg, etc.'
        item_list = [1, 2, 3]
        storage = CachedArray(typecode='I')
        storage.extend(item_list)
        assert storage.max == max(item_list)
        assert storage.min == min(item_list)
        assert storage.average == 2
        assert len(storage) == len(item_list)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_sequence_length_distrib']
    unittest.main()
