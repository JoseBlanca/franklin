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
from biolib.collections import SparseVector

class SparseVectorTests(unittest.TestCase):
    '''It tests the sparse vector collection'''
    @staticmethod
    def test_sparse_vector():
        'It tests the sparse vector collection'
        spv = SparseVector(nelements=100)
        spv[50] = 30
        assert spv[50] == 30

    @staticmethod
    def test_store_non_int():
        'sparse vectors can hold non int values'
        spv = SparseVector(nelements=100, store_non_int=True)
        spv[50] = [30]
        assert spv[50] == [30]


if __name__ == "__main__":
    unittest.main()
