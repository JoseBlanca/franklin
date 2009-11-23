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
from biolib.utils.misc_utils import (xml_itemize, _get_xml_tail,
                                     _get_xml_header, NamedTemporaryDir,
                                     FileIndex, split_long_sequences)
from biolib.utils.collections_ import FileCachedList
from biolib.seq.seqs import SeqWithQuality

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
        dir_name = temp_dir.name
        assert os.path.exists(dir_name) == True
        temp_dir.close()
        assert os.path.exists(dir_name) == False

        temp_dir = NamedTemporaryDir()
        dir_name = temp_dir.name
        fhand = open(os.path.join(dir_name, 'peio'), 'w')
        assert os.path.exists(fhand.name) == True
        assert os.path.exists(dir_name)   == True
        del(temp_dir)
        assert os.path.exists(dir_name) == False


class TestFileIndexer(unittest.TestCase):
    'It test the FileIndex class'

    @staticmethod
    def test_basic_file_index():
        'It test the file index class basic functionality'
        fhand = StringIO.StringIO('>key1\nhola\n>key2\ncaracola\n')
        index = FileIndex(fhand, item_start_patterns=['>'],
                          key_patterns=['>([^ \t\n]+)'])
        assert index['key1'] == '>key1\nhola\n'
        assert index['key2'] == '>key2\ncaracola\n'

    @staticmethod
    def test_with_item_types():
        'It test the file index class basic functionality'
        content = '>\n%key1%\ntype1\nhola\n>\n%key2%\ncaracola\ntype2\n'
        fhand = StringIO.StringIO(content)
        index = FileIndex(fhand, item_start_patterns=['>'],
                          key_patterns=['%([^%]+)%'],
                          type_patterns={'type1':['type1'], 'type2':['type2']})
        assert index['type1']['key1'] == '>\n%key1%\ntype1\nhola\n'
        assert index['type2']['key2'] == '>\n%key2%\ncaracola\ntype2\n'

    @staticmethod
    def test_item_types_from_file():
        'It test the file index class basic functionality'
        content = '>\n%key1%\ntype1\nhola\n>\n%key2%\ncaracola\ntype2\n'
        fhand = StringIO.StringIO(content)
        index = FileIndex(fhand, item_start_patterns=['>'],
                          key_patterns=['%([^%]+)%'],
                          type_patterns=['(type[0-9])'])
        assert index['type1']['key1'] == '>\n%key1%\ntype1\nhola\n'
        assert index['type2']['key2'] == '>\n%key2%\ncaracola\ntype2\n'



class SplitLongSequencestest(unittest.TestCase):
    'It tests sequence spliting functions'
    @staticmethod
    def test_split_long_sequences():
        '''It test the function that splits sequences of an iterator with long
        sequences into  smaller sequences'''
        seq = 'atatatatatg'
        qual = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        seq_rec = SeqWithQuality(seq=seq, qual=qual)
        seq_iter = iter([seq_rec])
        splited_seq_iter = split_long_sequences(seq_iter, 5)
        seq1 = splited_seq_iter.next()
        seq2 = splited_seq_iter.next()
        assert len(seq1) == 6
        assert len(seq2) == 5

        seq_iter = iter([seq_rec])
        splited_seq_iter = split_long_sequences(seq_iter, 3)
        seq1 = splited_seq_iter.next()
        seq2 = splited_seq_iter.next()
        seq3 = splited_seq_iter.next()
        assert len(seq1) == 4
        assert len(seq2) == 4
        assert len(seq3) == 3

class FileCachedListTest(unittest.TestCase):
    'It tests a list like class cached on a file'
    @staticmethod
    def test_filecachedlist():
        'It test the functionality of this list like class'
        #with ints
        clist = FileCachedList(type_=int)
        clist.append(0)
        clist.append(1)
        for index, item in enumerate(clist.items()):
            assert item == index
        #with floats
        clist = FileCachedList(type_=float)
        clist.append(0)
        clist.append(1)
        for index, item in enumerate(clist.items()):
            assert item == float(index)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
