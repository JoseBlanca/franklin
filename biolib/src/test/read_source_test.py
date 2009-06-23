'''
Created on 2009 eka 22

@author: peio
'''
import unittest
from biolib.read_source import (ReadSourceRegex, ReadSourceFile, 
                                get_read_strain, ReadSources)
from StringIO import StringIO

REGEX_LIST = [('\w\wLE\w*', 'comav_library1'),
              ('\w\wUH\w*', 'comav_library1')
                 ]
CLONE_READ_LIBRARY = '''read,clone,library
ESIE1234121,121313132clone, a
ESEE1234122,121313133clone,b
'''
class LibrarySourceRegexTest(unittest.TestCase):
    'It uses this test to LibrarySourceRegex'
    def test_basic_use(self):
        '''It test the basic use of the class ReadSourceRegex'''
        read_name = 'ESLE1234121'
        library_regex = ReadSourceRegex(REGEX_LIST)
        library_name = library_regex.get_library(read_name)
        assert  library_name == 'comav_library1'
        try:
            library_regex.get_clone(read_name)
            self.fail()
        except Exception:
            pass




class LibrarySourceFileTest(unittest.TestCase):
    'It uses this test to check LibrarySourceFile'
    @staticmethod
    def test_basic_use():
        '''It test the basic use of the class  ReadSourceFile'''
        read             = 'ESIE1234121'
        fhand_clone_read = StringIO(CLONE_READ_LIBRARY)
        read_source      = ReadSourceFile(fhand_clone_read)
        library          = read_source.get_library(read)
        clone            = read_source.get_clone(read)
        assert library == 'a'
        assert clone   == '121313132clone'
class ReadSourcesTest(unittest.TestCase):
    'It test ReadSources class'
    @staticmethod
    def test_basic_use():
        '''It test the basic use of the class ReadSources'''
        read             = 'ESLE1234121'
        read1            = 'ESEE1234122'
        fhand_clone_read = StringIO(CLONE_READ_LIBRARY)
        read_source      = ReadSourceFile(fhand_clone_read)
        library_regex    = ReadSourceRegex(REGEX_LIST)
        read_sources     = ReadSources([read_source, library_regex])
        library = read_sources.get_library(read)
        assert library == 'comav_library1'
        library = read_sources.get_library(read1)
        assert library == 'b'

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

class GetStrainTest(unittest.TestCase):
    'It test get_strain function'
    @staticmethod
    def test_basic_use():
        'It test get_strain function'
        library_name = 'a'
        fhand_library = StringIO(LIBRARY_FILE)
        strain = get_read_strain(library_name, [fhand_library])
        assert strain == 'Oregon-R'
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()