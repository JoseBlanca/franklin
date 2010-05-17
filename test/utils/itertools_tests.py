'''
Created on 11/05/2010

@author: peio
'''
import unittest
from franklin.utils.itertools_ import take_sample, make_cache, store

class TakeSampleTest(unittest.TestCase):
    'tests take sample test'
    @staticmethod
    def test_take_sample():
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

class MakeCacheTest(unittest.TestCase):
    'It test the cache for iterables'
    @staticmethod
    def test_make_cache():
        'It test the cache for iterables'
        item_list = [1,3, 4, 4.7, 'hola', [1, {'caracola':True}]]
        items = make_cache(iter(item_list))
        assert list(items) == item_list

class StoreTest(unittest.TestCase):
    'It test the store for iterables'
    @staticmethod
    def test_store():
        'It test the store for iterables'
        item_list = [1,3, 4, 4.7, 'hola', [1, {'caracola':True}]]
        storage = store()
        storage.extend(item_list)
        assert list(storage) == item_list

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
