'''
Created on 16/03/2010

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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
from tempfile import NamedTemporaryFile

from os.path import join

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.backbone.create_project import create_project
from franklin.backbone.analysis import BACKBONE_BASENAMES, BACKBONE_DIRECTORIES
from franklin.backbone.blast_runner import backbone_blast_runner, \
    guess_blastdb_kind

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
        query = '>seq1\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGCT'
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
        test_dir.close()

    def test_blast_seq_against_bad_db(self):
        'We can blast a seq file against a database'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = join(test_dir.name, project_name)

        #some query fasta file
        query = '>seq1\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGCT'
        query += 'CATACCCCTGCCGAACCGCTTTTGTCA|n'
        query_fhand = NamedTemporaryFile(mode='w')
        query_fhand.write(query)
        query_fhand.flush()

        #the blast db
        blast_db_fname = 'uni'
        blast_db = join(DATA_DIR, 'blast', blast_db_fname)

        blast_program = 'blastn'
        try:
            backbone_blast_runner(query_fpath=query_fhand.name,
                                  project_dir=project_dir,
                                  blast_program=blast_program,
                                  blast_db=blast_db)
            self.fail('RuntimeError expected')
        except RuntimeError:
            pass
        test_dir.close()

    @staticmethod
    def test_blast_seq_against_seq_db():
        'We can blast a seq file against a sequence file database'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        create_project(directory=test_dir.name,
                                       name=project_name)
        project_dir = join(test_dir.name, project_name)

        #some query fasta file
        query = '>seq1\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGCT'
        query += 'CATACCCCTGCCGAACCGCTTTTGTCA\n'
        query_fhand = NamedTemporaryFile(mode='w')
        query_fhand.write(query)
        query_fhand.flush()

        #the blast db
        blast = '@seq\nGATCGGCCTTCTTGCGCATCTCACGCGCTCCTGCGGCGGCCTGTAGGGCAGGC\n'
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
        assert '<Hit_accession>seq</Hit_accession>' in open(blast_fpath).read()

    @staticmethod
    def test_blastdb_seq_kind():
        'It test the blastdb kind'
        blastdb = join(DATA_DIR, 'blast', 'tomato_genome2')
        assert  guess_blastdb_kind(blastdb) == 'nucl'
        blastdb = '/srv/databases/blast/tair7_pep'
        assert  guess_blastdb_kind(blastdb) == 'prot'

if    __name__ == "__main__":
    #import sys;sys.argv = ['', 'BlastTest.test_blast_seq_against_seq_db']
    unittest.main()
