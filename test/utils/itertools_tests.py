'''
Created on 11/05/2010

@author: peio
'''
import unittest
from franklin.utils.itertools_ import (take_sample, make_cache, store, classify,
                                       ungroup)
import itertools

class TakeSampleTest(unittest.TestCase):
    'tests take sample test'
    @staticmethod
    def test_take_sample():
        'tests take sample test'
        #basic test
        items = iter(range(100))
        sample = list(take_sample(items, 10))
        assert len(sample) == 10
        assert sample != [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        n_items = [1000, 100, 100000, 100000, 100000]
        sample_sizes = [990, 8, 100, 100, 200]

        for n_item, sample_size in zip(n_items, sample_sizes):
            repeats = 0
            while repeats < 10:
                repeats += 1
                iterator= iter(range(n_item))
                a = take_sample(iterator, sample_size)
                assert sample_size ==  len(list(a))

    def test_tee_sample(self):
        'It tests that tee and sample behave ok together'
        items = iter(range(1000))
        sample = take_sample(items, 50)
        sample1, sample2 = itertools.tee(sample)

class MakeCacheTest(unittest.TestCase):
    'It tests the cache for iterables'
    @staticmethod
    def test_make_cache():
        'It test the cache for iterables'
        item_list = [1, 3, 4, 4.7, 'hola', [1, {'caracola':True}]]
        items = make_cache(iter(item_list))
        assert list(items) == item_list

class StoreTest(unittest.TestCase):
    'It tests the store for iterables'
    @staticmethod
    def test_store():
        'It test the store for iterables'
        item_list = [1, 3, 4, 4.7, 'hola', [1, {'caracola':True}]]
        storage = store()
        storage.extend(item_list)
        assert list(storage) == item_list

class ClassifierTest(unittest.TestCase):
    'It tests the classifier function'
    @staticmethod
    def test_classifier():
        'It tests the classifier function'
        item_list = [1, 0.5, 2, .75, 'a']
        classification = classify(item_list, lambda x: type(x).__name__)
        assert list(classification['int']) == [1, 2]
        assert list(classification['float']) == [0.5, 0.75]
        assert list(classification['str']) == ['a']

class UngroupTest(unittest.TestCase):
    'It tests the ungorup generator'
    @staticmethod
    def test_ungroup():
        'It tests the ungroup generator'
        items = [{'a': 1, 'b':2}, {'a':2}, {'c':3}]
        items = list(ungroup(items, lambda x: x.items()))
        assert items == [('a', 1), ('b', 2), ('a', 2), ('c', 3)]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
