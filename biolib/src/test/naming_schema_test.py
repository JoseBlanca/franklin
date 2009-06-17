'''It test that we can modify the names/accessions in different kind of files
like ace, caf or fasta'''

import unittest
import os.path
from StringIO import StringIO
import biolib
from biolib.naming_schema import (change_names_in_files,
                                  create_names_for_contigs,
                                  create_naming_database, DbNamingSchema,
                                  FileNamingSchema,
                                  add_project_to_naming_database)
from tempfile import NamedTemporaryFile
import sqlalchemy

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
'''),
'library':('''
format-version:1
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
    properties: SO:222, SO:3456
''')
}

class NamingSchema(object):
    'A mock Naming schema class'
    #pylint: disable-msg=W0613
    #pylint: disable-msg=R0903
    def __init__(self, **kwargs):
        '''It won't matter what you pass it will always return 001, 002, etc.'''
        self._current = '000'
    def get_uniquename(self, name, kind):
        '''An alias for get_next_name.'''
        number =  str(int(self._current) + 1).rjust(3, '0')
        self._current = number
        return kind + number
    def commit_last_name(self):
        'The pega'
        pass

class DbNamingSchemaTest(unittest.TestCase):
    'We test the behaviour of a naming schema stored in a naming database'
    
    @staticmethod
    def test_basic_behaviour():
        'It tests that we can get names'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(engine)
        add_project_to_naming_database(engine, name='my_project', code='my',
                                       description='a test project')
        naming = DbNamingSchema(engine, project='my_project')
        assert naming.get_uniquename(kind='EST') == 'myES000001'
        assert naming.get_uniquename(kind='EST') == 'myES000002'
        assert naming.get_uniquename(kind='transcribed_cluster') == 'myTC000001'
        naming.commit()

    @staticmethod
    def test_name_persistance():
        'The names are stored between db commits'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(engine)
        add_project_to_naming_database(engine, name='my_project', code='my',
                                       description='a test project')
        naming = DbNamingSchema(engine, project='my_project')
        assert naming.get_uniquename(kind='EST') == 'myES000001'
        naming.commit()
        naming = DbNamingSchema(engine, project='my_project')
        assert naming.get_uniquename(kind='EST') == 'myES000002'
        
    @staticmethod
    def test_rollback():
        "If we don't commit we lose the changes"
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(engine)
        add_project_to_naming_database(engine, name='my_project', code='my',
                                       description='a test project')
        naming = DbNamingSchema(engine, project='my_project')
        assert naming.get_uniquename(kind='EST') == 'myES000001'

        naming = DbNamingSchema(engine, project='my_project')
        assert naming.get_uniquename(kind='EST') == 'myES000001'

class FileNamingSchemaTest(unittest.TestCase):
    'We test the behaviour of a naming schema cached in a file'
    
    def test_basic_behaviour(self):
        'It tests that we can get names using a FileNamingSchema'
        fhand = NamedTemporaryFile()
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(engine)
        add_project_to_naming_database(engine, name='my_project', code='my',
                                       description='a test project')
        naming = DbNamingSchema(engine, project='my_project')
        naming = FileNamingSchema(fhand, naming)
        assert naming.get_uniquename(kind='EST', name='hola') == 'myES000001'
        assert naming.get_uniquename(name='hola') == 'myES000001'
        assert naming.get_uniquename(kind='EST', name='caracol') == 'myES000002'
        naming.commit()
        fhand.seek(0)
        naming = FileNamingSchema(fhand)
        assert naming.get_uniquename(name='hola') == 'myES000001'
        assert naming.get_uniquename(name='caracol') == 'myES000002'
        try:
            assert naming.get_uniquename(kind='EST', name='pascual')
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass

    def test_roolback(self):
        'If the changes are not commited they are ignored'
        fhand = NamedTemporaryFile()
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(engine)
        add_project_to_naming_database(engine, name='my_project', code='my',
                                       description='a test project')
        naming = DbNamingSchema(engine, project='my_project')
        naming = FileNamingSchema(fhand, naming)
        assert naming.get_uniquename(kind='EST', name='hola') == 'myES000001'
        naming.commit()
        assert naming.get_uniquename(kind='EST', name='caracol') == 'myES000002'

        naming = FileNamingSchema(fhand)
        try:
            assert naming.get_uniquename(name='caracol')
            self.fail()
            #pylint: disable-msg=W0704
        except ValueError:
            pass

#class ChangeNameTest(unittest.TestCase):
#    '''It test that we can modify the names/accs in different kind of files'''
#    @staticmethod
#    def test_fasta():
#        'It test that we can change the name in the fasta files.'
#        fhand_in = StringIO(EXAMPLES['fasta'][0])
#        fhand_out = StringIO('')
#        naming = NamingSchema('aproject', 'feature_kind',
#                              'a db connection')
#        naming = CachedNamingSchema(naming)
#        change_names_in_files(fhand_in, fhand_out, naming, 'fasta')
#        assert fhand_out.getvalue() == EXAMPLES['fasta'][1]
#
#    @staticmethod
#    def test_contig_naming():
#        'It tests that we can create the names for a caf file.'
#        #caf naming
#        contig_path = os.path.join(DATA_DIR, 'example.caf')
#        contig_fh = open(contig_path, 'rt')
#        est_naming = NamingSchema('aproject', 'est', 'a db connection')
#        contig_naming = NamingSchema('aproject', 'contig', 'a db connection')
#        names = create_names_for_contigs(contig_fh, 'caf', {'read':est_naming,
#                                                      'contig':contig_naming})
#        assert names['Contig1']       == '001'
#        assert names['22ak65e10.s1t'] == '001'
#        assert names['22ak93c2.s1t']  == '021'
#
#        #ace naming
#        contig_path = os.path.join(DATA_DIR, 'example.ace')
#        contig_fh = open(contig_path, 'rt')
#        est_naming = NamingSchema('aproject', 'est', 'a db connection')
#        contig_naming = NamingSchema('aproject', 'contig', 'a db connection')
#        names = create_names_for_contigs(contig_fh, 'ace', {'read':est_naming,
#                                                      'contig':contig_naming})
#        assert names['Contig2']          == '001'
#        assert names['eucalyptus0111']   == '001'
#        assert names['eucalyptus57514']  == '034'
#
#        #caf and ace complete name change
#        for kind in ('ace', ('caf')):
#            fhand_in = StringIO(EXAMPLES[kind][0])
#            fhand_out = StringIO('')
#            est_naming = NamingSchema('aproject', 'est', 'a db connection')
#            contig_naming = NamingSchema('aproject', 'contig',
#                                         'a db connection')
#            names = create_names_for_contigs(fhand_in, kind,
#                                             {'read':est_naming,
#                                              'contig':contig_naming})
#            change_names_in_files(fhand_in, fhand_out, names, kind)
#            assert fhand_out.getvalue() == EXAMPLES[kind][1]

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testiprscan_parse']
    unittest.main()
    
