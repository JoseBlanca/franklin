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
import StringIO
import biolib
from biolib.fpc import gff_parser, FPCMap, write_fpc_gff, fpcgff2_parser

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

    @staticmethod
    def test_gff():
        'It tests the fpc gff3 writting'
        fhand = tempfile.TemporaryFile(suffix='.gff')
        fpc_fname = os.path.join(DATA_DIR, 'fpc_test.fpc')
        fpc = FPCMap(open(fpc_fname))
        write_fpc_gff(fpc, fhand=fhand)

    @staticmethod
    def test_fpcgff2_parsing():
        'It tests fpcgff2_parser'
        pfcgff2file = '''
Chr0\tassembly\tChromosome\t1\t1091372440\t.\t.\t.\tSequence "Chr0"; Name "Chr0"
Chr0\tFPC\tcontig\t20357139\t20983827\t.\t.\t.\tcontig "ctg19"; Name "ctg19"
Chr0\tFPC\tBAC\t20549651\t20934675\t.\t.\t.\tBAC "Cm27_D06"; Name "Cm27_D06"; Contig_hit "19"
Chr0\tFPC\tBAC\t20406291\t20840467\t.\t.\t.\tBAC "Cm14_J14"; Name "Cm14_J14"; Marker_hit "A_21-C11 0 0"; Contig_hit "19"
Chr0\tFPC\tBAC\t20434963\t20967443\t.\t.\t.\tBAC "Cm32_M24"; Name "Cm32_M24"; Marker_hit "A_21-C11 0 0"; Contig_hit "19"
Chr0\tFPC\tBAC\t20496403\t20938771\t.\t.\t.\tBAC "Cm03_D18"; Name "Cm03_D18"; Contig_hit "19"
Chr0\tFPC\tBAC\t20455443\t20983827\t.\t.\t.\tBAC "Cm06_I11"; Name "Cm06_I11"; Contig_hit "19"
Chr0\tFPC\tBAC\t20451347\t20856851\t.\t.\t.\tBAC "Cm32_P11"; Name "Cm32_P11"; Contig_hit "19"
Chr0\tFPC\tBAC\t20357139\t20963347\t.\t.\t.\tBAC "Cm27_G15"; Name "Cm27_G15"; Marker_hit "A_21-C11 0 0"; Contig_hit "19"
Chr0\tFPC\tmarker\t20701203\t20701203\t.\t.\t.\tmarker "A_21-C11"; Name "A_21-C11"; Contig_hit "ctg19 - 3" (Cm14_J14 Cm32_M24 Cm27_G15)
'''

        pfcgff2fhand = StringIO.StringIO(pfcgff2file)

        features = list(fpcgff2_parser(pfcgff2fhand))



        assert features[0]['name'] == 'Chr0'
        assert features[1]['name'] == "ctg19"
        assert features[1]['parents'] == ['Chr0']
        assert features[2]['name'] == "Cm27_D06"
        assert features[2]['parents'] == ['ctg19']
        assert features[9]['name']  == "A_21-C11"
        assert "ctg19" in features[9]['parents']
        assert "Cm14_J14" in features[9]['parents']
        assert "Cm27_G15" in features[9]['parents']

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestFPC.test_fpc']
    unittest.main()