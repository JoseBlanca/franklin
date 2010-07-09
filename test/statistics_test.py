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
import os, random, tempfile

from franklin.statistics import (CachedArray, histogram, create_distribution,
                                 draw_boxplot)
import franklin

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

class DistributionTest(unittest.TestCase):
    'It tests the create distribution function'
    def test_basic_distribution(self):
        'It tests the distribution'
        summary_fhand = StringIO()
        distrib_fhand = StringIO()

        numbers = CachedArray(typecode='I')
        numbers.extend([1, 2, 3, 4, 5, 6, 7, 8, 9, 101, 2, 3, 4, 5, 6, 7, 8, 9])

        create_distribution(numbers, distrib_fhand=distrib_fhand,
                            summary_fhand=summary_fhand)
        result = '''minimum: 1
maximum: 101
average: 10.56
variance: 486.91
sum: 190
items: 18'''
        assert result in summary_fhand.getvalue()

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

class BoxPlotTest(unittest.TestCase):
    'It checks the boxplot drawing'

    def test_boxplot(self):
        'It checks the boxplot drawing'
        #some data
        mu = 10
        sigma = 3
        num_values = 1000
        lists = []
        for index in range(5):
            values = [random.normalvariate(mu, sigma)
                                                for index_ in range(num_values)]
            lists.append(values)
            values = [random.uniform(mu+sigma, mu-sigma)
                                                for index_ in range(num_values)]
            lists.append(values)
        plot_fhand = tempfile.NamedTemporaryFile()
        stats_fhand = StringIO()
        draw_boxplot(lists, xlabel='distributions', ylabel='distrib',
                     title='boxplot', fhand=open(plot_fhand.name, 'w'),
                     stats_fhand=stats_fhand)
        assert 'PNG' in  open(plot_fhand.name).read(10)
        plot_fhand.close()

        result = stats_fhand.getvalue()
        assert '09' in result
        assert 'median' in result

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
    #import sys;sys.argv = ['', 'BoxPlotTest.test_boxplot']
    unittest.main()
