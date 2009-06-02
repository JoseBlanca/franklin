'''It test that we can modify the names/accessions in different kind of files
like ace, caf or fasta'''

import unittest
import os.path
from StringIO import StringIO
import biolib
from biolib.naming_schema import (change_names_in_files, CachedNamingSchema,
                                  create_names_for_contigs)

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
Assembled_from EST1

DNA : contig1
ACGTACTGTAG

BaseQuality : contig1

DNA : EST1
ATCGATCTGTAC

BaseQuality : EST1
Sequence : EST1
Is_read
''', '''
Sequence : contig1
Is_contig
Assembled_from EST1

DNA : contig1
ACGTACTGTAG

BaseQuality : contig1

DNA : EST1
ATCGATCTGTAC

BaseQuality : EST1
Sequence : EST1
Is_read
'''
)

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

class ChangeNameTest(unittest.TestCase):
    '''It test that we can modify the names/accs in different kind of files'''
    @staticmethod
    def test_fasta():
        'It test that we can change the name in the fasta files.'
        #for kind in EXAMPLES:
        for kind in ['fasta']:
            fhand_in = StringIO(EXAMPLES[kind][0])
            fhand_out = StringIO('')
            naming = NamingSchema('aproject', 'feature_kind',
                                  'a db connection')
            naming = CachedNamingSchema(naming)
            change_names_in_files(fhand_in, fhand_out, naming, 'fasta')
            assert fhand_out.getvalue() == EXAMPLES[kind][1]

    @staticmethod
    def test_contig_naming():
        'It tests that we can create the names for a caf file.'
        #caf
        contig_path = os.path.join(DATA_DIR, 'example.caf')
        contig_fh = open(contig_path, 'rt')
        est_naming = NamingSchema('aproject', 'est', 'a db connection')
        contig_naming = NamingSchema('aproject', 'contig', 'a db connection')
        names = create_names_for_contigs(contig_fh, 'caf', {'read':est_naming,
                                                      'contig':contig_naming})
        assert names['Contig1']       == '001'
        assert names['22ak65e10.s1t'] == '001'
        assert names['22ak93c2.s1t']  == '021'

        #ace
        contig_path = os.path.join(DATA_DIR, 'example.ace')
        contig_fh = open(contig_path, 'rt')
        est_naming = NamingSchema('aproject', 'est', 'a db connection')
        contig_naming = NamingSchema('aproject', 'contig', 'a db connection')
        names = create_names_for_contigs(contig_fh, 'ace', {'read':est_naming,
                                                      'contig':contig_naming})
        print names
        assert names['Contig2']          == '001'
        assert names['eucalyptus0111']   == '001'
        assert names['eucalyptus57514']  == '034'

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
    
