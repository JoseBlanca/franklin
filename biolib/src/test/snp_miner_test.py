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
from biolib.seqvar.snp_miner import (SNPMINER_MAP_DEF,
                                      create_snp_miner_database,
                                      add_snv_to_db, add_reference_annot)
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
    def test_add_snp_to_db():
        '''It test the simple case - add snp'''
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)

        ref1 = 'ref1'
        lib1_alleles = {'library': 'some_libary',
                        'annotations': {},
                        'alleles': [{'allele':'A', 'reads':3, 'kind':INVARIANT,
                                     'qualities':[30, 30, 30]},
                                    {'allele':'T', 'reads':3, 'kind':SNP,
                                     'qualities':[30, 30, 30]},]
                        }
        snv = Snv(name='snp1', reference=ref1, location=3,
                      lib_alleles=[lib1_alleles])

        add_snv_to_db(engine, snv, library='lib1')
        snp_miner = DbMap(engine, SNPMINER_MAP_DEF)
        ref_inst = snp_miner.select_one('reference', {'name':'ref1'})
        assert ref_inst.name == 'ref1'
        snp_inst  = snp_miner.select_one('seqvar', {'name':'snp1'})
        assert snp_inst.name ==  'snp1'
        assert snp_inst.location.reference.name == 'ref1'

    @staticmethod
    def xtest_add_reference_annot():
        '''It test the simple case - add contig_annot'''
        annotation_dict = {'name': 'ref1', 'type':'blast',
                            'value':'GO:0022'}

        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)

        ref1 = SeqWithQuality(seq='atatatat', name='ref1')
        seq_var = SeqVariation(name='snp1', alleles= [{'allele':'A', 'reads':3,
                                                       'kind':INVARIANT,
                                                       'quality':[30, 30, 30]},
                                                       {'allele':'T', 'reads':3,
                                                        'kind':SNP,
                                                       'quality':[30, 30, 30]}],
                                reference=ref1, location=3)

        add_seqvar_to_db(engine, seq_var, library='lib1')
        add_reference_annot(engine, [annotation_dict])
        snp_miner = DbMap(engine, SNPMINER_MAP_DEF)
        ref_inst = snp_miner.select_one('referenceprop', {'type':'blast'})
        ref_inst2 = snp_miner.select_one('reference', {'name':'ref1'})
        assert ref_inst.reference_id == ref_inst2.reference_id

if __name__ == "__main__":
    unittest.main()


