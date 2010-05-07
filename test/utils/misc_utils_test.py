'''
Created on 2009 mai 22

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
import StringIO
from franklin.utils.misc_utils import (xml_itemize, _get_xml_tail,
                                       _get_xml_header, NamedTemporaryDir,
                                       VersionedPath, get_num_threads)
from franklin.utils.collections_ import FileCachedList

class Minor_utilities_test(unittest.TestCase):
    'Test form minor utilities'
    @staticmethod
    def test_get_num_threads():
        'tests get_num_threads'
        threads = 3
        assert get_num_threads(threads) == threads

        threads = False
        assert get_num_threads(threads) == 1


class XMLTest(unittest.TestCase):
    '''It tests the xml utils'''

    @staticmethod
    def test_xml_itemize():
        '''It tests xml itemize'''
        string = '<h><t><c></c><c></c></t></h>'
        xml = StringIO.StringIO(string)
        cont = 0
        for result in  xml_itemize(xml, 'c'):
            assert result == '<h><t><c></c></t></h>'
            cont += 1
        assert cont == 2

    @staticmethod
    def test_xml_itemize_by_chunks():
        '''It tests xml itemize but with more than one item'''
        string = '<h><t><c>1</c><c>2</c><c>3</c><c>4</c><c>5</c></t></h>'
        xml = StringIO.StringIO(string)
        xmls = list(xml_itemize(xml, 'c', num_items=2))
        assert xmls[0] == '<h><t><c>1</c><c>2</c></t></h>'
        assert xmls[1] == '<h><t><c>3</c><c>4</c></t></h>'
        assert xmls[2] == '<h><t><c>5</c></t></h>'
        assert len(xmls) == 3

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
        assert os.path.exists(dir_name) == True
        del(temp_dir)
        assert os.path.exists(dir_name) == False

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

class VersionedPathTest(unittest.TestCase):
    'It tests the versioned path class'
    def test_basic_functionality(self):
        'VersionedPath basic functionality'
        tempdir = NamedTemporaryDir()
        tempdir_name = tempdir.name
        fnames = ['hola.txt', 'hola.0.txt', 'hola.1.txt',
                  'adios.txt', 'adios.1.txt',
                  'foo.txt']
        [open(os.path.join(tempdir.name, fname), 'w') for fname in fnames]

        path_str = os.path.join(tempdir.name, 'hola.txt')
        path = VersionedPath(path_str)
        assert str(path) == path_str

        path_str_0 = os.path.join(tempdir.name, 'hola.0.txt')
        path = VersionedPath(path_str_0)
        assert str(path) == path_str
        assert path.basename == 'hola'
        assert path.directory == tempdir.name
        assert path.extension == 'txt'

        assert path.last_version == VersionedPath(os.path.join(tempdir_name,
                                                               'hola.1.txt'))
        assert path.next_version == VersionedPath(os.path.join(tempdir_name,
                                                               'hola.2.txt'))

        fpaths = [os.path.join(tempdir_name, fname) for fname in ('hola.1.txt',
                                                                  'adios.1.txt',
                                                                  'foo.txt')]
        expected_paths = set(fpaths)
        versioned_paths = set(path.list_fpaths_versioned())
        assert versioned_paths == expected_paths

        path = VersionedPath(tempdir_name)
        versioned_paths = set(path.list_fpaths_versioned())
        assert versioned_paths == expected_paths

        tempdir.close()

        #in an empty dir
        path = VersionedPath('hola.txt')
        assert path.last_version.endswith('hola.txt')
        assert path.next_version.endswith('hola.0.txt')

        path = VersionedPath('pl_illumina.sm_rp_75_59_uc82.sfastq')
        assert path.last_version.endswith('pl_illumina.sm_rp_75_59_uc82.sfastq')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
