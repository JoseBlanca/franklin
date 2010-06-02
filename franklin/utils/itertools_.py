'''
Created on 04/09/2009

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# project is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

#import pysparse
from __future__  import division
import tempfile

import itertools, random
from array import array
import cPickle as pickle

from tempfile import TemporaryFile

def list_consecutive_pairs_iter(items):
    'Given a list it yields all consecutive pairs (1,2, 2,3, 3,4)'
    first_time = True
    first_item, second_item = None, None
    for item in items:
        if first_time:
            first_item = item
            first_time = False
            continue
        second_item = item
        yield first_item, second_item
        first_item, second_item = second_item, first_item

def take_sample(iterator, sample_size, num_items_in=None):
    '''This function takes a sample of size sample size from the given iterator.

    Optionaly the function offers the possibility to give the num of items in
    the iterator to improve velocity
    '''

    if num_items_in is None:
        iterator, iterator2 = itertools.tee(iterator, 2)
        num_items_in = 0
        for item in iterator2:
            num_items_in += 1

    if num_items_in < sample_size:
        return iterator
    else:
        return _take_sample(iterator, sample_size, num_items_in)

def _take_sample(iterator, sample_size, num_items_in):
    '''This function takes a sample of size sample size from the given iterator.

    This is the function that do the job
    '''
    rate = sample_size / num_items_in

    reserved_rate = rate
    if reserved_rate * num_items_in < 10:
        reserved_rate = 10 / num_items_in
    reserve_length = int(num_items_in * reserved_rate)

    counter = 0
    reserved_items = []
    reserve_full = False
    for item in iterator:
        random_num = random.random()
        if random_num < rate:
            yield item
            counter += 1
            if counter >= sample_size:
                break
        elif not reserve_full and random_num - rate < reserved_rate:
            reserved_items.append(item)
            if len(reserved_items) >= reserve_length:
                reserve_full = True
    items_left = sample_size - counter
    for item in random.sample(reserved_items, items_left):
        yield item

def make_cache(iterator):
    'Given an iterator it makes a new iterator that feeds from a cache'

    cache_fhand = TemporaryFile(suffix='.cache')
    #we save the old iterator in the cache
    for item in iterator:
        pickle.dump(item, cache_fhand)

    #now we create the generator that feeds from the cache
    cache_fhand.seek(0)
    while True:
        try:
            yield pickle.load(cache_fhand)
        except EOFError:
            break

    cache_fhand.close() #we remove the temporary cache

class store():
    '''It stores pickable elements for future use.

    Similar to a queue but once you ask for the first element you can add no
    more.
    '''
    def __init__(self):
        'It inits the store'
        self._cache_fhand = tempfile.TemporaryFile()
        self._lock = False

    def append(self, item):
        'It adds one item to the store'
        if not self._lock:
            pickle.dump(item, self._cache_fhand)
        else:
            raise RuntimeError('The store is locked')

    def extend(self, items):
        'It adds the items to the store'
        if not self._lock:
            for item in items:
                self.append(item)
        else:
            raise RuntimeError('The store is locked')

    def next(self):
        'It returns one item'
        if not self._lock:
            self._lock = True
            self._cache_fhand.seek(0)
        try:
            return pickle.load(self._cache_fhand)
        except EOFError:
            raise StopIteration

    def __iter__(self):
        'Part of the iterator protocol'
        return self

def classify(items, classifier):
    '''Given an iterator and a classifier function it returns several iterators

    It returns a dict. They keys will be the values returned by the classifier
    when is applied to each item.
    The items should be pickable and the values returned by the classifier
    hashable.
    '''
    classifications = {}
    for item in items:
        kind = classifier(item)
        if kind not in classifications:
            classifications[kind] = store()
        classifications[kind].append(item)
    return classifications

def ungroup(items, ungrouper):
    '''Similar to map, but the ungrouper returns an iterable

    All items returned by the ungrouper will be added to the generator.
    '''
    for item in items:
        subitems = ungrouper(item)
        for subitem in subitems:
            yield subitem

MIN_FREE_MEMORY_PERCENT = 10
MEMORY_CHECK_CYCLES = 10000

def _check_free_memory(percent=None):
    'It checks that we still have more free memory than the given percent'
    if not percent:
        percent = MIN_FREE_MEMORY_PERCENT
    fhand = open('/proc/meminfo', 'r')
    mem_total = int(fhand.readline().split()[1])
    mem_free = int(fhand.readline().split()[1])
    free_percent = mem_free / mem_total * 100
    if percent > free_percent:
        raise RuntimeError('Scarce memory')

class CachedArray(object):
    '''It stores numbers for future use.

    It can work in memory or write it contents to a file.
    '''
    def __init__(self, typecode):
        'The init '
        self._typecode = typecode
        self._array = array(typecode)
        self._cache_fhand = None
        self._max = None
        self._min = None
        self._len = 0
        self._sum = 0
        self._check_in_cycles = MEMORY_CHECK_CYCLES

    def extend(self, items):
        'It adds all items to the store'
        for item in items:
            self.append(item)

    def append(self, item):
        'It appends one item to the store'
        #some statistics
        if self._max is None:
            self._max = item
            self._min = item
        else:
            if self._max < item:
                self._max = item
            if self._min > item:
                self._min = item
        self._sum += item
        self._len += 1

        #are we running low on memory?
        self._check_in_cycles -= 1
        if self._check_in_cycles <= 0 and self._cache_fhand is not None:
            try:
                _check_free_memory()
            except RuntimeError:
                #we are running low on memory, so we store in disk
                self.to_disk()
            self._check_in_cycles = MEMORY_CHECK_CYCLES

        if self._array is not None:
            self._array.append(item)
        else:
            self._cache_fhand.write(str(item) + '\n')

    def _generate_file_items(self):
        'It yields all items from the file cache'
        #let's be sure that all items are in the disk
        if self._cache_fhand:
            self._cache_fhand.flush()

        casts = {'h': int, 'H':int, 'i':int, 'I':int, 'L':int, 'l':int,
                 'f':float }
        cast = casts[self._typecode]
        for line in open(self._cache_fhand.name):
            yield cast(line)

    def _get_items(self):
        'A generator that yields all items'
        if self._array is not None:
            return iter(self._array)
        else:
            return self._generate_file_items()

    def __iter__(self):
        'Part of the iterator protocol'
        return self._get_items()

    def to_disk(self):
        'It saves all items to the disk'
        if self._cache_fhand is not None:
            return  #we are already in disk
        self._cache_fhand = tempfile.NamedTemporaryFile()
        #we store all item in the file
        for item in self._array:
            self._cache_fhand.write(str(item) + '\n')
        self._cache_fhand.flush()
        #we remove the array storage
        del self._array
        self._array = None

    def _get_min(self):
        'It returns the minimum'
        return self._min
    min = property(_get_min)

    def _get_max(self):
        'It returns the maximum'
        return self._max
    max = property(_get_max)

    def _get_sum(self):
        'It returns the sum'
        return self._sum
    sum = property(_get_sum)

    def _get_average(self):
        'It returns the maximum'
        return self._sum / self._len
    average = property(_get_average)

    def __len__(self):
        'It returns the number of items'
        return self._len

    def _get_variance(self):
        'It returns the variance'
        sum_ = 0
        count = 0
        mean = self.average
        for number in self:
            sum_ += (number - mean) ** 2
            count += 1
        if not count:
            return None
        else:
            return (sum_ / float(count))
    variance = property(_get_variance)
