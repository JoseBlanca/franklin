'''
Created on 2009 eka 19

@author: peio
'''
from django.views import static

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
from biolib.seqvar.snp_miner import (SNPMINER_MAP_DEF, SnvDb,
                                     create_snp_miner_database)

from biolib.seqvar.seqvariation import Snv, INVARIANT, SNP
from biolib.seqs import SeqWithQuality
from biolib.db.db_utils import DbMap
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
        kind      = SNP
        snv_sql = snp_miner.get_snv_sql(reference, location, kind)
        assert snv_sql.reference.name == 'ref1'
        assert snv_sql.location == 2

        snv_sql = snp_miner.get_snv_sql(reference, location, kind)
        assert  snv_sql.reference.name == 'ref1'

    @staticmethod
    def test_add_alleles_per_library():
        'It tests add_alleles_per_library method'
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
        kind      = INVARIANT
        snv_sql = snv_miner.get_snv_sql(reference, location, kind)

        snv_miner.add_alleles_per_library(snv_sql, library_info1)
        allele_sql = snv_miner.select_one('LibrarySnvAlleles', {'allele':'A'})
        assert allele_sql.kind == INVARIANT
        assert allele_sql.reads == 3
        assert eval(allele_sql.qualities) == [20, 30, 30]

        #the snv has to be changed by now
        assert snv_sql.kind == SNP

        library_sql = snv_miner.select_one('Library',
                                         {'accession':library_info1['library']})
        assert library_sql.accession == library_info1['library']

    @staticmethod
    def test_add_annotation_per_library():
        'It tests add_alleles_per_library method'
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        snv_miner = SnvDb(engine)

        ref1 = 'ref1'
        library_info1 = {'library': 'some_libary',
                        'annotations': {'hola': 'caracola'}}

        # add snv_sql
        reference = ref1
        location  = 2
        kind      = INVARIANT
        snv_sql = snv_miner.get_snv_sql(reference, location, kind)

        snv_miner.add_annotations_per_library(snv_sql, library_info1)
        annot = snv_miner.select_one('LibrarySnvAnnots', {'kind':'hola'})
        assert annot.value == 'caracola'

if __name__ == "__main__":
    unittest.main()


