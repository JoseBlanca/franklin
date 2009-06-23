'''
Created on 2009 eka 22

@author: peio
'''
import unittest
from StringIO import StringIO 
from biolib.file_parsers import library_parser

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