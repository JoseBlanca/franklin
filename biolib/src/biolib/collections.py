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
class SparseVector(object):
    'A sparse vector collection to hold nearly empty vectors.'
    def __init__(self, nelements, size_hint=None, store_non_int=False):
        '''It initializes the vector.

        We require the number of elements that the vector should hold.
        size_hint is the number of estimated elements that the vector should
        hold at the end. If given expensive memory reallocation will be saved.
        '''
        if size_hint:
            self._spv = pysparse.spmatrix.ll_mat(nelements, 1, size_hint)
        else:
            self._spv = pysparse.spmatrix.ll_mat(nelements, 1)
        if store_non_int:
            self._cargo = []
            self._cargo.append(None)   #we dont' use the 0 as cargo index
        else:
            self._cargo = None

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
