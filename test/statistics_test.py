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
                                 draw_boxplot, draw_histogram,
                                 draw_stacked_columns, IntsStats)
import franklin

TEST_DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'data')

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
        result = '''Statistics for histogram
-------------------------
minimum: 1
maximum: 101
average: 10.5556
variance: 486.9136
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
        numpy_distrib = ([2, 2, 4, 2, 3], [0, 2, 4, 6, 8, 10])
        our_distrib = histogram(numbers, bins=bins)
        for num1, num2 in zip(numpy_distrib[0], our_distrib[0]):
            assert num1 == num2
        for num1, num2 in zip(numpy_distrib[1], our_distrib[1]):
            assert num1 == num2
            assert isinstance(num2, int)

        #with floats
        numbers = [0.0, 1, 2, 3, 4, 4, 4, 5, 6, 7, 8, 9, 10]
        bins    = 5
        numpy_distrib = ([2, 2, 4, 2, 3], [0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
        our_distrib = histogram(numbers, bins=bins)

        for num1, num2 in zip(numpy_distrib[0], our_distrib[0]):
            assert num1 == num2
        for num1, num2 in zip(numpy_distrib[1], our_distrib[1]):
            assert num1 == num2
            assert isinstance(num2, float)

        #now we draw the histogram
        #draw_histogram(our_distrib[0], our_distrib[1])
        png_file = tempfile.NamedTemporaryFile(suffix='.png')
        #draw_histogram(our_distrib[0], our_distrib[1], fhand=png_file)
        #assert 'PNG' in open(png_file.name).read(4)
        svg_file = tempfile.NamedTemporaryFile(suffix='.svg')
        draw_histogram(our_distrib[0], our_distrib[1], fhand=svg_file)
        svg_file.flush()
        assert 'xml' in open(svg_file.name).read(5).lower()

    @staticmethod
    def test_remove_outliers():
        'It tests the histogram removing the outliers'
        numbers = [random.gauss(5, 1) for index in range(100)]
        numbers.append(5000)
        bins = 10
        distrib = histogram(numbers, bins=bins, remove_outliers=True)
        assert int(distrib[1][-1]) < 20

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
        plot_fhand = tempfile.NamedTemporaryFile(suffix='.svg')
        stats_fhand = StringIO()
        draw_boxplot(lists, xlabel='distributions', ylabel='distrib',
                     title='boxplot', fhand=plot_fhand,
                     stats_fhand=stats_fhand)
        assert 'xml' in  open(plot_fhand.name).read(10)

        result = stats_fhand.getvalue()
        assert '09' in result
        assert 'median' in result
        plot_fhand.close()


        #the not to draw all boxes
        plot_fhand = tempfile.NamedTemporaryFile(suffix='.svg')
        draw_boxplot(lists, xlabel='distributions', ylabel='distrib',
                     title='boxplot', fhand=plot_fhand,
                     stats_fhand=stats_fhand, max_plotted_boxes=5)
        assert 'xml' in  open(plot_fhand.name).read(10)
        plot_fhand.close()

class DrawFreqHistogram(unittest.TestCase):
    'It test the freq histogram'

    @staticmethod
    def test_stacked_columns():
        'It test the freq histogram'
        values = {
              'A': [0.66666666666666663, 0.33333333333333331, 0.0, 0.0, 0, 0.0],
              'C': [0.33333333333333331, 0.66666666666666663, 0.0, 0.0, 0, 0.0],
              'T': [0.0, 0.0, 1.0, 0.33333333333333331, 0, 1.0],
              'G': [0.0, 0.0, 0.0, 0.66666666666666663, 0, 0.0]}

        colors = {'A':'g', 'C':'b', 'T':'r', 'G':'k'}
        fhand = tempfile.NamedTemporaryFile(suffix='.svg')
        draw_stacked_columns(values, title='test', fhand=fhand, colors=colors)
        fhand.flush()
#        from franklin.utils.cmd_utils import call
#        call(('eog', fhand.name))
        svg = open(fhand.name, 'r').read()
        assert '<!-- Created with matplotlib (http://matplotlib' in svg


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

    @staticmethod
    def test_sample():
        'It tests the random sample'
        storage = CachedArray(typecode='I')
        storage.extend([0] * 10000)
        storage.extend([1] * 10000)
        storage.extend([2] * 10000)
        count = {}
        for item in storage.sample:
            try:
                count[item] += 1
            except KeyError:
                count[item] = 1
        len_sample = len(storage.sample)
        assert len_sample == storage.sample_length

        assert count[1] / (len_sample * 1.0) - 1/3.0 < 0.05

        assert storage.get_sample_percentiles([5, 95]) == [0, 2]


class IntsStatsTest(unittest.TestCase):
    'It test the extensible array class'
    @staticmethod
    def create_test_array():
        d = {'9':'5', '10':'288', '11':'002556688', '12':'00012355555',
             '13':'0000013555688', '14':'00002555558',
             '15':'0000000000355555555557', '16':'000045', '17':'000055',
             '18':'0005', '19':'00005', '21':'5'}
        ext_array = IntsStats()
        for key, values in d.items():
            for num in values:
                ext_array.append(int(key+num))
        return ext_array

    @staticmethod
    def test_array():
        'Create an extensible array'
        ext_array = IntsStats(init_len=5)
        ext_array.append(6)
        ext_array.append(2)
        assert  ext_array.min == 2
        assert  ext_array.max == 6
        ext_array.append(200)
        assert ext_array.max == 200

        input_ = (3, 5, 7, 7, 38)
        ext_array = IntsStats(input_)
        assert ext_array.median == 7

    def test_distribution(self):
        'It tests the histogram function'
#        ints_array = IntsStats(1)
#        num_integers = 1000000
#        max_value    = 100000
#        for i in xrange(num_integers):
#            ints_array.append(random.randint(0, max_value))

#        (distrib, bin_edges) = ints_array.histogram(25)
#
#        (distrib, bin_edges) = ints_array.histogram(25, range_=None)
#
#        (distrib, bin_edges) = ints_array.histogram(25, range_=(None, 1000))
#
#        (distrib, bin_edges) = ints_array.histogram(25, range_=(1000, None))
#
#        (distrib, bin_edges) = ints_array.histogram(25, range_=(1000, 2000))

        ints_array = self.create_test_array()
        distrib = ints_array.calculate_distribution(bins=10,
                                                        remove_outliers=5)

        assert distrib['distrib'] == [7L, 13L, 7L, 10L, 7L, 22L, 6L, 4L, 5L, 5L]
        assert distrib['bin_edges'] == [110, 118, 126, 134, 142, 150, 158, 166,
                                        174, 182, 190]

        general_stats = '''Statistics for histogram
-------------------------
minimum: 95
maximum: 215
average: 145.1522
variance: 557.4334
sum: 13354
items: 92'''
        result_fhand = tempfile.NamedTemporaryFile()
        ints_array.write_general_stats(result_fhand)
        assert general_stats in open(result_fhand.name).read()

        ints_array = IntsStats([0, 0, 1, 3])

        assert [2, 1, 1]  == ints_array.calculate_distribution(bins=3)['distrib']

    def test_stats_functs(self):
        'It test the statistical functions of the class'
        ext_array = IntsStats()
        ext_array.append(3)
        ext_array.append(5)
        ext_array.append(7)
        ext_array.append(7)
        ext_array.append(38)
        assert ext_array.median == 7

        ext_array = IntsStats()
        ext_array.append(3)
        ext_array.append(5)
        ext_array.append(7)
        ext_array.append(7)
        assert ext_array.median == 6

        ext_array = self.create_test_array()
        assert ext_array.median == 145
        assert round(ext_array.average, 2) == 145.15

        assert ext_array.sum == 13354
        assert ext_array.count == 92
        assert round(ext_array.variance, 2) == 557.43


if __name__ == "__main__":
#    import sys;sys.argv = ['', 'IntsStatsTest.test_distribution']
    unittest.main()
