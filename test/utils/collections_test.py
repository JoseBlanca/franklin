'''
Created on 11/05/2010

@author: peio
'''
import unittest
from franklin.utils.collections_ import FileCachedList, take_sample

class FileCachedListTest(unittest.TestCase):
    'It tests a list like class cached on a file'
    @staticmethod
    def test_filecachedlist():
        'It test the functionality of this list like class'
        #with ints
        clist = FileCachedList(type_=int)
        clist.append(0)
        clist.append(1)
        for index, item in enumerate(clist.items()):
            assert item == index
        #with floats
        clist = FileCachedList(type_=float)
        clist.append(0)
        clist.append(1)
        for index, item in enumerate(clist.items()):
            assert item == float(index)

class TakeSampleTest(unittest.TestCase):
    'tests take sample test'
    @staticmethod
    def test_take_sample_test():
        'tests take sample test'
        n_items = [1000, 100, 100000, 100000, 100000]
        sample_sizes = [990, 8, 100, 100, 200]

        for n_item, sample_size in zip(n_items, sample_sizes):
            repeats = 0
            while repeats < 10:
                repeats += 1
                iterator= iter(range(n_item))
                a = take_sample(iterator, sample_size)
                assert sample_size ==  len(list(a))









if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
