'''
Created on 2009 eka 22

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

import unittest
from StringIO import StringIO
from franklin.gmod.file_parsers import library_parser

LIBRARY_FILE = '''format-version:1
library_definition
    name: a
    type: library type:genomic
    organism:Cucumis melo
    cvterms: SO:0001, SO:0002
    properties: property type:strain:Oregon-R, property type:stage:adult male

library_definition
    name:b
    type: library type:genomic
    organism: Cucumis melo
    cvterms:SO:0003, SO:0004
    properties: property type:strain:a_fly, property type:stage:pupa
'''

class FileParserTest(unittest.TestCase):
    '''Tests file parsers '''
    @staticmethod
    def test_library_parser():
        '''It test library parser'''
        fhand = StringIO(LIBRARY_FILE)
        cont = 0
        for library in library_parser(fhand):
            cont += 1
        assert  cont == 2


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
