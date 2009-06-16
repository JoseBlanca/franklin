'''It test that we can modify the names/accessions in different kind of files
like ace, caf or fasta'''

import unittest
import os.path
from StringIO import StringIO
import biolib
from biolib.naming_schema import (change_names_in_files, CachedNamingSchema,
                                  create_names_for_contigs,
                                  generate_naming_file)

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

EXAMPLES = {'fasta':('''
>hola
ATCGTAGTCAGTTCAGTCTATGCTAGT
>carcola
ATGCGTAGCTAGTCGTAGTCTAGTCAT''','''
>001
ATCGTAGTCAGTTCAGTCTATGCTAGT
>002
ATGCGTAGCTAGTCGTAGTCTAGTCAT'''),
'caf': ('''
Sequence : contig1
Is_contig
Assembled_from EST1 1 293 39 331

DNA : contig1
ACGTACTGTAG

BaseQuality : contig1

DNA : EST1
ATCGATCTGTAC

BaseQuality : EST1
Sequence : EST1
Is_read''',
'''
Sequence : 001
Is_contig
Assembled_from 001 1 293 39 331

DNA : 001
ACGTACTGTAG

BaseQuality : 001

DNA : 001
ATCGATCTGTAC

BaseQuality : 001
Sequence : 001
Is_read'''),
'ace':('''
CO Contig2 213 6 1 U
TCGTGGTGGTTCATCGACATT
        
BQ
96 97 95 95 97 93 96

AF eucalyptus0111 U 1

RD eucalyptus0111 173 0 0
TCGTGGTGGTTCATCGAC

QA 1 173 1 173
''', '''
CO 001 213 6 1 U
TCGTGGTGGTTCATCGACATT
        
BQ
96 97 95 95 97 93 96

AF 001 U 1

RD 001 173 0 0
TCGTGGTGGTTCATCGAC

QA 1 173 1 173
''')
}

class NamingSchema(object):
    'A mock Naming schema class'
    #pylint: disable-msg=W0613
    #pylint: disable-msg=R0903
    def __init__(self, project, feature_kind, db_connection):
        '''It won't matter what you pass it will always return 001, 002, etc.'''
        self._current = '000'
    def get_next_name(self):
        'It returns the next name'
        number =  str(int(self._current) + 1).rjust(3, '0')
        self._current = number
        return number
    def __getitem__(self, name):
        '''An alias for get_next_name.'''
        return self.get_next_name()
    def commit_last_name(self):
        'The pega'
        pass

class ChangeNameTest(unittest.TestCase):
    '''It test that we can modify the names/accs in different kind of files'''
    @staticmethod
    def test_fasta():
        'It test that we can change the name in the fasta files.'
        fhand_in = StringIO(EXAMPLES['fasta'][0])
        fhand_out = StringIO('')
        naming = NamingSchema('aproject', 'feature_kind',
                              'a db connection')
        naming = CachedNamingSchema(naming)
        change_names_in_files(fhand_in, fhand_out, naming, 'fasta')
        assert fhand_out.getvalue() == EXAMPLES['fasta'][1]

    @staticmethod
    def test_contig_naming():
        'It tests that we can create the names for a caf file.'
        #caf naming
        contig_path = os.path.join(DATA_DIR, 'example.caf')
        contig_fh = open(contig_path, 'rt')
        est_naming = NamingSchema('aproject', 'est', 'a db connection')
        contig_naming = NamingSchema('aproject', 'contig', 'a db connection')
        names = create_names_for_contigs(contig_fh, 'caf', {'read':est_naming,
                                                      'contig':contig_naming})
        assert names['Contig1']       == '001'
        assert names['22ak65e10.s1t'] == '001'
        assert names['22ak93c2.s1t']  == '021'

        #ace naming
        contig_path = os.path.join(DATA_DIR, 'example.ace')
        contig_fh = open(contig_path, 'rt')
        est_naming = NamingSchema('aproject', 'est', 'a db connection')
        contig_naming = NamingSchema('aproject', 'contig', 'a db connection')
        names = create_names_for_contigs(contig_fh, 'ace', {'read':est_naming,
                                                      'contig':contig_naming})
        assert names['Contig2']          == '001'
        assert names['eucalyptus0111']   == '001'
        assert names['eucalyptus57514']  == '034'

        #caf and ace complete name change
        for kind in ('ace', ('caf')):
            fhand_in = StringIO(EXAMPLES[kind][0])
            fhand_out = StringIO('')
            est_naming = NamingSchema('aproject', 'est', 'a db connection')
            contig_naming = NamingSchema('aproject', 'contig',
                                         'a db connection')
            names = create_names_for_contigs(fhand_in, kind,
                                             {'read':est_naming,
                                              'contig':contig_naming})
            change_names_in_files(fhand_in, fhand_out, names, kind)
            assert fhand_out.getvalue() == EXAMPLES[kind][1]
            
class GeneratingNamingFileTest(unittest.TestCase):
    ''' It tests the generate_naming_file function'''
    @staticmethod
    def test_generate_naming_file():
        ''' It tests the generate_naming_file function'''
        
        naming_schema  = NamingSchema('a_prohject', 'lib', 'a connection') 
        naming_schemas = {}
        naming_schemas['library_names'] = naming_schema
        fpath = os.path.join(DATA_DIR, 'library.txt')
        fhand = open(fpath, 'r')
        generate_naming_file(fhand, naming_schemas, 'library' )
    @staticmethod
    def test_generate_naming_with_2():
        '''It test the naming change with two naming schemas'''
        contig_name_schema = NamingSchema('a_prohject', 'lib', 'a connection')
        read_naming_schema = NamingSchema('a_prohject', 'lib', 'a connection')
        naming_schemas = {}
        naming_schemas['read_names']   = read_naming_schema
        naming_schemas['contig_names'] = contig_name_schema
        fpath = os.path.join(DATA_DIR, 'example.caf')
        fhand = open(fpath, 'r')
        generate_naming_file(fhand, naming_schemas, 'caf' )
    def test_eror_generate_naming_file(self):
        ''' It tests error detection of the generate_naming_file function'''
        
        naming_schema  = NamingSchema('a_prohject', 'lib', 'a connection') 
        naming_schemas = {}
        naming_schemas['library_n'] = naming_schema
        fpath = os.path.join(DATA_DIR, 'library.txt')
        fhand = open(fpath, 'r')
        try:
            generate_naming_file(fhand, naming_schemas, 'library' )
            self.fail()
        except ValueError:
            pass
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
    
