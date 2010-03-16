'''
Created on 16/03/2010

@author: jose
'''

import unittest, os
from tempfile import NamedTemporaryFile

from os.path import join

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.backbone.create_project import create_project
from franklin.backbone.analysis import BACKBONE_BASENAMES, BACKBONE_DIRECTORIES
from franklin.backbone.blast_runner import backbone_blast_runner

def _get_basename(fpath):
    'It returns the base name without path and extension'
    return os.path.splitext(os.path.basename(fpath))[0]

class BlastTest(unittest.TestCase):
    'It test the blast runner'

    @staticmethod
    def test_blast_seq_against_db():
        'We can blast a seq file against a database'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = join(test_dir.name, project_name)

        #some query fasta file
        query  = '>seq1\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGCT'
        query += 'CATACCCCTGCCGAACCGCTTTTGTCA|n'
        query_fhand = NamedTemporaryFile(mode='w')
        query_fhand.write(query)
        query_fhand.flush()

        #the blast db
        blast_db_fname = 'univec+'
        blast_db = join(DATA_DIR, 'blast', blast_db_fname)

        blast_program = 'blastn'
        backbone_blast_runner(query_fpath=query_fhand.name,
                              project_dir=project_dir,
                              blast_program=blast_program,
                              blast_db=blast_db)

        #is the blast ok?
        blast_fpath = join(project_dir,
                           BACKBONE_DIRECTORIES['blast_dir'],
                           _get_basename(query_fhand.name),
                           blast_db_fname,
                           '%s.%s.xml' % (BACKBONE_BASENAMES['blast_basename'],
                                          blast_program))
        assert '<Hit_def>vec1</Hit_def>' in open(blast_fpath).read()

    @staticmethod
    def test_blast_seq_against_seq_db():
        'We can blast a seq file against a sequence file database'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = join(test_dir.name, project_name)

        #some query fasta file
        query  = '>seq1\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGCT'
        query += 'CATACCCCTGCCGAACCGCTTTTGTCA\n'
        query_fhand = NamedTemporaryFile(mode='w')
        query_fhand.write(query)
        query_fhand.flush()

        #the blast db
        blast  = '@seq\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGC\n'
        blast += '+\n'
        blast += '11111111111111111111111111111111111111111111111111111\n'
        blast_db_fhand = NamedTemporaryFile(mode='w', suffix='.sfastq')
        blast_db_fhand.write(blast)
        blast_db_fhand.flush()

        blast_program = 'blastn'
        backbone_blast_runner(query_fpath=query_fhand.name,
                              project_dir=project_dir,
                              blast_program=blast_program,
                              blast_db_seq=blast_db_fhand.name)

        #is the blast ok?
        blast_fpath = join(project_dir,
                           BACKBONE_DIRECTORIES['blast_dir'],
                           _get_basename(query_fhand.name),
                           _get_basename(blast_db_fhand.name),
                           '%s.%s.xml' % (BACKBONE_BASENAMES['blast_basename'],
                                          blast_program))
        assert '<Hit_def>seq</Hit_def>' in open(blast_fpath).read()

if    __name__    ==    "__main__":
    #import sys;sys.argv = ['', 'BlastTest.test_blast_seq_against_seq_db']
    unittest.main()
