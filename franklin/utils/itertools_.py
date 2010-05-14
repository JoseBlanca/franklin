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

import itertools, random

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
    rate = num_items_in / sample_size

    reserved_rate = 0.01 if 0.005 * num_items_in > 10 else 10 / num_items_in

    counter = 0
    reserved_items = []
    for item in iterator:
        random_num = random.random()
        if random_num < rate :
            yield item
            counter += 1
            if counter >= sample_size:
                break
        elif random_num < reserved_rate:
            reserved_items.append(item)
    else:
        for reserved_item in reserved_items:
            if counter <= sample_size:
                yield reserved_item
                counter += 1
            else:
                break
