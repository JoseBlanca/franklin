'''
Created on 2009 mai 22

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

import unittest, os
import StringIO
from biolib.biolib_utils import (xml_itemize, _get_xml_tail, _get_xml_header,
                                NamedTemporaryDir, guess_seq_file_format)
class XMLTest(unittest.TestCase):
    '''It tests the xml utils'''

    @staticmethod
    def test_xml_itemize():
        '''It tests xml itemize '''
        string = '<h><t><c></c><c></c></t></h>'
        xml = StringIO.StringIO(string)
        cont = 0
        for result in  xml_itemize(xml, 'c'):
            assert result == '<h><t><c></c></t></h>'
            cont += 1
        assert cont == 2
    def test_no_good_xml_start_end(self):
        '''Tests if the raise an error with a bad xml file. from begining to
        end '''
        xml = StringIO.StringIO('<header><conten></content></header>')
        self.failUnlessRaises(ValueError, _get_xml_header, xml, 'content')
    def test_no_good_xml_end_start(self):
        '''Tests if the raise an error with a bad xml file. From end to start'''
        xml = StringIO.StringIO('<header><content><content></header>')
        self.failUnlessRaises(ValueError, _get_xml_tail , xml, 'content')

class NamedTemporariDirTest(unittest.TestCase):
    'It test temporay named dir'
    @staticmethod
    def test_simple_named_temporary_dir():
        'It test temporay named dir'
        temp_dir = NamedTemporaryDir()
        dir_name = temp_dir.name()
        assert os.path.exists(dir_name) == True
        temp_dir.close()
        assert os.path.exists(dir_name) == False

        temp_dir = NamedTemporaryDir()
        dir_name = temp_dir.name()
        fhand = open(os.path.join(dir_name, 'peio'), 'w')
        assert os.path.exists(fhand.name) == True
        assert os.path.exists(dir_name)   == True
        del(temp_dir)
        assert os.path.exists(dir_name) == False

class GuessFormatSeqFileTest(unittest.TestCase):
    'It tests that we can guess the format of a sequence file'
    @staticmethod
    def test_guess_format():
        'It test that we can guess the format for the sequence files'
        fhand = StringIO.StringIO('>fasta\n')
        assert guess_seq_file_format(fhand) == 'fasta'

        fhand = StringIO.StringIO('LOCUS AX0809\n')
        assert guess_seq_file_format(fhand) == 'genbank'

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
