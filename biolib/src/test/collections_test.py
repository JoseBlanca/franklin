'''Tests for the collections module'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

import unittest
from biolib.collections_ import SparseVector
from biolib.rbtree import RBDict

class SparseVectorTests(unittest.TestCase):
    '''It tests the sparse vector collection'''
    @staticmethod
    def test_sparse_vector():
        'It tests the sparse vector collection'
        spv = RBDict()
        spv[50] = 30
        spv[40] = 30
        assert [(40, 30), (50, 30)] == list((spv.items()))

    @staticmethod
    def test_store_non_int():
        'sparse vectors can hold non int values'
        spv = SparseVector(nelements=4, store_non_int=True)
        spv[0] = 30
        spv[2] = [30]
        assert spv[2] == [30]
        assert spv[0] == 30
        result = [30, None, [30], None]
        for index, item in enumerate(spv):
            assert item == result[index]
        spv[1] = 'a'
        assert list(spv.non_empty_items()) == [(0, 30), (1, 'a'), (2, [30])]
        assert list(spv.non_empty_values()) == [30, 'a', [30]]

if __name__ == "__main__":
    unittest.main()
