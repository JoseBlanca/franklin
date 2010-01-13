'''It test that we can modify the names/accessions in different kind of files
like ace, caf or fasta'''

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

import unittest
import os.path
from StringIO import StringIO
import biolib
from biolib.db.naming import (change_names_in_files,
                              create_naming_database, DbNamingSchema,
                              FileNamingSchema,
                              add_project_to_naming_database,
                              project_in_database)
from tempfile import NamedTemporaryFile
import sqlalchemy
from biolib.utils.misc_utils import DATA_DIR

EXAMPLES = {'fasta':('''
>hola
ATCGTAGTCAGTTCAGTCTATGCTAGT
>carcola
ATGCGTAGCTAGTCGTAGTCTAGTCAT''','''
>myES000001
ATCGTAGTCAGTTCAGTCTATGCTAGT
>myES000002
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
'''),
'fastq':('''
@MEES000003
AGATGGATGTTCCACAGGCTTTGATATGGTCCTCGACTTACTACGCAAATCGATGCCAGATCCTCCAATCCAGATATATT
CAGCAATAATCTTGTCGGTATATGGAGTGACGT
+
EEEEEEEEEA@@@@DDD@@444@@IFEC:;:;;AEEEDDEEEEEEA===EEEEEEEEEEEEEEDDDDDEEEEEEEEEEEE
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEDA
@MEES000004
CTGTCTCTGCGAAGAAGAGCAAGGCATCCAGTAGTTCGGAGGAAGACTCATCTGAAGATGATTCTGATTCTGATGAGGAA
CCAAAGGGAAAGTTATTGCAGCAA
+
FFFFFFFFFFFFFFFFFFIIIIIIIIIIIIIIIIIIIIIIIIIIIFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
FFBBBF;;>>78@FEEEEFFFFFF
''')
}

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

    @staticmethod
    def test_project_in_database():
        'It test the ability to know if the databse is already added'
        project = 'test'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(engine)
        assert not project_in_database(engine, project)
        add_project_to_naming_database(engine, name=project, code='my',
                                       description='a test project')
        assert project_in_database(engine, project)
        assert not project_in_database(engine, 'test2')

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

class ChangeNameTest(unittest.TestCase):
    '''It test that we can modify the names/accs in different kind of files'''
    @staticmethod
    def test_fasta():
        'It test that we can change the name in the fasta files.'
        fhand_in  = StringIO(EXAMPLES['fasta'][0])
        fhand_out = StringIO('')
        engine    = sqlalchemy.create_engine('sqlite:///:memory:')

        create_naming_database(engine)
        add_project_to_naming_database(engine, name='test_project', code='my',
                                       description='a test project')
        naming    = DbNamingSchema(engine, 'test_project')

        change_names_in_files(fhand_in, fhand_out, naming, 'fasta', 'EST')
        assert fhand_out.getvalue() == EXAMPLES['fasta'][1]


    @staticmethod
    def test_fastq():
        'It test that we can change the name in the fasta files.'
        fhand_in  = open(os.path.join(DATA_DIR, 'solexa.fastq'))
        fhand_out = StringIO('')
        engine    = sqlalchemy.create_engine('sqlite:///:memory:')

        create_naming_database(engine)
        add_project_to_naming_database(engine, name='test_project', code='my',
                                       description='a test project')
        naming    = DbNamingSchema(engine, 'test_project')

        change_names_in_files(fhand_in, fhand_out, naming, 'fastq', 'EST')
        assert fhand_out.getvalue()[:8] == '@myES000'

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

