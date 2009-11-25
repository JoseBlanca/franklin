'''
Created on 23/09/2009

@author: jose
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

import unittest, os, tempfile

import biolib
from biolib.gmod.fpc import FPCMap, write_fpc_gff

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

class TestFPC(unittest.TestCase):
    'It tests the fpc functionality'

    @staticmethod
    def test_fpc():
        'It tests the fpc parsing'
        fpc_fname = os.path.join(DATA_DIR, 'fpc_test.fpc')
        fpc = FPCMap(open(fpc_fname))
        assert fpc.name == 'demo'
        assert fpc.version == '8.5.1'
        assert len(fpc.markers) == 2
        assert len(fpc.clones) == 4
        assert len(fpc.contigs) == 2
        #print fpc.contigs

    @staticmethod
    def test_gff():
        'It tests the fpc gff3 writting'
        fhand = tempfile.TemporaryFile(suffix='.gff')
        fpc_fname = os.path.join(DATA_DIR, 'fpc_test.fpc')
        fpc = FPCMap(open(fpc_fname))
        write_fpc_gff(fpc, fhand=fhand)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestFPC.test_fpc']
    unittest.main()