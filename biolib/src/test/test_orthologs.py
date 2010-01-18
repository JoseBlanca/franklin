'''
Created on 15/01/2010

@author: peio
'''
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

import os, unittest
from biolib.utils.misc_utils import DATA_DIR
from biolib.orthologs import get_orthologs

class OrthologsTests(unittest.TestCase):
    'It test basic use of the orthologs functions'
    @staticmethod
    def test_get_orthologs():
        'It tests the get_orthologs functions'
        blast_file  = open(os.path.join(DATA_DIR, 'melon_tair.xml'))
        blast_file2 = open(os.path.join(DATA_DIR, 'tair_melon.xml'))

        orthologs = get_orthologs(blast_file, blast_file2)
        #print orthologs.next()
        assert orthologs.next() == ('melon1', 'tair1')
        assert orthologs.next() == ('melon2', 'tair2')



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
