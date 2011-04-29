'''
Created on 2009 mai 21

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

import unittest, os
from franklin.iprscan import xml_iprscan_parser_iter
from franklin.utils.misc_utils import TEST_DATA_DIR

class IprsacnTest(unittest.TestCase):
    '''It test all functions in iprscan parser '''

    @staticmethod
    def test_iprscan_parse():
        '''Tests iprscan parser '''
        fname = os.path.join(TEST_DATA_DIR, 'iprscan.xml')
        fhand = open(fname, 'r')
        protein_ids = ["PRRB_MYCTU", "Q9RHD9_PSEAE", "RS16_ECOLI"]
        cont = 0
        for protein_result in xml_iprscan_parser_iter(fhand):
            assert  protein_result['id'] in protein_ids
            cont += 1
        assert cont == 3



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
