'''
Created on 26/10/2009

@author: jose
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

import unittest
from StringIO import StringIO

from franklin.gff import write_gff

class GffOutTest(unittest.TestCase):
    'It tests the gff writer'

    @staticmethod
    def test_simple_output():
        'We can write a simple gff3 file'
        feat1 = {'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 'id':    'gene00001',
                 'name':  'EDEN'
                 }
        feats = ['##sequence-region ctg123 1 1497228', feat1]
        result = '''##gff-version 3
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tID=gene00001;Name=EDEN\n'''
        outh = StringIO()
        write_gff(feats, outh)
        assert result == outh.getvalue()

        feat1 = {'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 'parents': ['p1', 'p2']}
        feats = [feat1]
        result = '''##gff-version 3
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tParent=p1,p2\n'''
        outh = StringIO()
        write_gff(feats, outh)
        assert result == outh.getvalue()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()