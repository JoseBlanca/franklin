'''
Created on 2009 eka 19

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

import unittest
from biolib.snv.snv_miner import SnvDb, create_snp_miner_database

from biolib.snv.snv import Snv, INVARIANT, SNP
import sqlalchemy

#from biolib.db.naming import (create_naming_database,
#                                  FileNamingSchema)

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


class SnpMinerTest(unittest.TestCase):
    '''It test the snp_miner package '''

    @staticmethod
    def test_get_snv_sql():
        'It test get_svn_sql function'

        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        snp_miner = SnvDb(engine)

        ref1 = 'ref1'

        reference = ref1
        location  = 2
        snv_sql = snp_miner.get_snv_sql(reference, location)
        assert snv_sql.reference.name == 'ref1'
        assert snv_sql.location == 2

        snv_sql = snp_miner.get_snv_sql(reference, location)
        assert  snv_sql.reference.name == 'ref1'

    @staticmethod
    def test_add_alleles_per_library():
        'It tests create_alleles_per_library method'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        snv_miner = SnvDb(engine)

        ref1 = 'ref1'
        library_info1 = {'library': 'some_libary',
                        'annotations': {'hola': 'caracola'},
                        'alleles': [{'allele':'A', 'reads':3, 'kind':INVARIANT,
                                     'qualities':[20, 30, 30]},
                                    {'allele':'T', 'reads':3, 'kind':SNP,
                                     'qualities':[30, 30, 30]},]}

        # add snv_sql
        reference = ref1
        location  = 2
        snv_sql = snv_miner.get_snv_sql(reference, location)
        snv_miner.create_alleles_per_library(snv_sql, library_info1)
        allele_sql = snv_miner.select_one('LibrarySnvAlleles', {'allele':'A'})
        assert allele_sql.kind == INVARIANT
        assert allele_sql.reads == 3
        assert eval(allele_sql.qualities) == [20, 30, 30]

        library_sql = snv_miner.select_one('Library',
                                         {'accession':library_info1['library']})
        assert library_sql.accession == library_info1['library']

    @staticmethod
    def test_add_annotation_per_library():
        'It tests create_alleles_per_library method'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        snv_miner = SnvDb(engine)

        ref1 = 'ref1'
        library_info1 = {'library': 'some_libary',
                        'annotations': {'hola': 'caracola'}}

        # add snv_sql
        reference = ref1
        location  = 2
        snv_sql = snv_miner.get_snv_sql(reference, location)

        snv_miner.create_annotations_per_library(snv_sql, library_info1)
        annot = snv_miner.select_one('LibrarySnvAnnots', {'kind':'hola'})
        assert annot.value == 'caracola'

    def test_add_snv_to_db(self):
        'It tests that we can add a snv to the database'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        snv_miner = SnvDb(engine)

        ref1 = 'ref1'
        library_info1 = {'library': 'some_library1',
                        'annotations': {'hola1': 'caracola1'},
                        'alleles': [{'allele':'A', 'reads':3, 'kind':INVARIANT,
                                     'qualities':[20, 30, 30]},
                                    {'allele':'T', 'reads':3, 'kind':SNP,
                                     'qualities':[30, 30, 30]},]}
        snv = Snv(reference=ref1, location=3,
                  per_lib_info=[library_info1])
        # add an snv
        snv_miner.create_snv(snv)
        selected_svn = snv_miner.select_one('Snv', attributes={'location':3})
        assert selected_svn.location == 3
        snv_miner.commit()
        # try to add again the same snv
        try:
            snv_miner.create_snv(snv)
            self.fail()
        except Exception:
            snv_miner.rollback()

        library_info2 = {'library': 'some_library2',
                        'annotations': {'hola': 'caracola'},
                        'alleles': [{'allele':'A', 'reads':3, 'kind':INVARIANT,
                                     'qualities':[20, 30, 30]},
                                    {'allele':'T', 'reads':3, 'kind':SNP,
                                     'qualities':[30, 30, 30]},]}
        snv = Snv(reference=ref1, location=3,
                  per_lib_info=[library_info2])

        # Trying to add another library to the same snv
        snv_miner.create_snv(snv)
        selected_snv = snv_miner.select_one('Snv', attributes={'location':3})
        lib_snvs = snv_miner.select('LibrarySnv',
                                    attributes={'snv': selected_snv})
        for lib_snv in lib_snvs:
            assert lib_snv.library.accession in ['some_library1',
                                                 'some_library2']
        assert selected_svn.location == 3

        # test the svn selection
        new_snv = snv_miner.select_snv(snv.reference, snv.location)

        assert new_snv.per_lib_info[0]['annotations'][0]['value'] == 'caracola1'
        assert new_snv.per_lib_info[1]['annotations'][0]['value'] == 'caracola'
        assert len(new_snv.per_lib_info) == 2

        #test snvs selections
        new_snv = list(snv_miner.select_snvs())[0]
        assert new_snv.per_lib_info[0]['annotations'][0]['value'] == 'caracola1'
        assert new_snv.per_lib_info[1]['annotations'][0]['value'] == 'caracola'
        assert len(new_snv.per_lib_info) == 2


if __name__ == "__main__":
    unittest.main()


