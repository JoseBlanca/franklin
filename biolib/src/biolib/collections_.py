'''
Created on 04/09/2009

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# project is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import pysparse
import tempfile

class SparseVector(object):
    'A sparse vector collection to hold nearly empty vectors.'
    def __init__(self, nelements, size_hint=None, store_non_int=False):
        '''It initializes the vector.

        We require the number of elements that the vector should hold.
        size_hint is the number of estimated elements that the vector should
        hold at the end. If given expensive memory reallocation will be saved.
        '''
        self._size = nelements
        if size_hint:
            self._spv = pysparse.spmatrix.ll_mat(nelements, 1, size_hint)
        else:
            self._spv = pysparse.spmatrix.ll_mat(nelements, 1)
        if store_non_int:
            self._cargo = []
            self._cargo.append(None)   #we dont' use the 0 as cargo index
        else:
            self._cargo = None
        self._position = 0

    def __setitem__(self, index, value):
        'It sets one item in the vector'
        cargo = self._cargo
        if cargo is None:
            self._spv[index, 0] = value
        else:
            #was an item already stored in that position?
            cargo_index = self._spv[index, 0]
            if cargo_index:
                self._cargo[cargo_index] = value
            else:
                self._cargo.append(value)
                cargo_index = len(self._cargo) - 1
                self._spv[index, 0] = cargo_index

    def __getitem__(self, index):
        'It returns one element'
        cargo = self._cargo
        if cargo is None:
            return self._spv[index, 0]
        else:
            return cargo[int(self._spv[index, 0])]
    def __iter__(self):
        'The iterable protocol'
        self._position = 0
        return self

    def next(self):
        'It returns the next item'
        position = self._position
        if position < self._size:
            value = self[self._position]
            self._position += 1
            return value
        else:
            raise StopIteration

    def __len__(self):
        'It returns the number of elements'
        return self._size

    def __str__(self):
        'It prints the vector'
        toprint = '['
        for item in self:
            toprint += str(item)
            toprint += ', '
        toprint += ']'
        return toprint

    def non_empty_items(self):
        'It return a iterator with non empty items'
        keys = self._spv.keys()[0]
        last_item = len(keys) - 1
        for index, key in enumerate(keys):
            if self._cargo:
                yield key, self._cargo[int(self._spv[key, 0])]
            else:
                yield key, self._spv[key, 0]
            if index == last_item:
                raise StopIteration

    def non_empty_values(self):
        'It return a iterator with non empty items'
        for item in self.non_empty_items():
            yield item[1]

class FileCachedList(object):
    '''A list cached in a file.

    The main aim is to store really big lists in files with an iterator
    interface.
    '''
    def __init__(self, type_):
        '''It creates a FileCachedList.

        It requires the type of objects that will be stored.
        '''
        self._type = type_
        #the file to store the data
        self._wfhand = tempfile.NamedTemporaryFile()
        #the file to read the data
        self._rfhand = None

    def append(self, item):
        'It writes one element at the file end'
        self._wfhand.write('%s\n' % str(item))

    def extend(self, items):
        'It adds a bunch of items from a list or from a FileCachedList'
        if 'items' in dir(items):
            items = items.items()
        for item in items:
            self.append(item)

    def items(self):
        'It yields all items'
        self._wfhand.flush()
        rfhand = open(self._wfhand.name)
        for line in rfhand:
            if len(line) == 0:
                raise StopIteration
            yield self._type(line.strip())
