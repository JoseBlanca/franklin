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
from biolib.db.snp_miner import (SNPMINER_MAP_DEF, create_snp_miner_database,
                                 add_contig_to_db, add_snp_to_db,
                                 add_contig_annot)
from biolib.db.db_utils import DbMap
import os, biolib, sqlalchemy
from StringIO import StringIO
from biolib.read_source import  ReadSourceFile

#from biolib.db.naming import (create_naming_database,
#                                  FileNamingSchema)

DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

SNP_FILE='''
format-version:1
snp
    name:snp1
    contig:Contig1
    start:3
    end:3
    kind:snp
    alleles: A:read1,read2 ; T:read3, read4
    annotations: pik:0.2
snp
    name:snp2
    contig:Contig1
    start:3
    end:3
    kind:snp
    alleles: A:read1,read2 ; T:read3, read4
    annotations: pik:0.2
'''
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
CLONE_READ_LIBRARY = '''read,clone,library
read1,121313132clone, a
read2,121313133clone,b
read3,121313134clone, a
read4,121313135clone,b
'''

class SnpMinerTest(unittest.TestCase):
    '''It test the snp_miner package '''
    @staticmethod
    def test_add_contig_to_db():
        '''It test the simple case - add contig'''
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        add_contig_to_db(engine, fhand)
        snp_miner = DbMap(engine, SNPMINER_MAP_DEF)
        snp_inst = snp_miner.select_one('contig', {'name':'Contig1'})
        assert snp_inst.name == 'Contig1'

    @staticmethod
    def test_add_snp_to_db():
        '''It test the simple case - add snp'''
        fhand_snp = StringIO(SNP_FILE)
        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        add_contig_to_db(engine, fhand)
        fhand_read_source = StringIO(CLONE_READ_LIBRARY)
        fhand_library = StringIO(LIBRARY_FILE)
        read_source = ReadSourceFile(fhand_read_source)
        add_snp_to_db(engine, fhand_snp, read_source, [fhand_library])
        snp_miner = DbMap(engine, SNPMINER_MAP_DEF)
        snp_inst  = snp_miner.select_one('snp', {'name':'snp1'})
        assert snp_inst.name ==  'snp1'
        assert snp_inst.contig.name == 'Contig1'
        snp_inst  = snp_miner.select_one('snpprop', {'snp_id':{'name':'snp1'},
                                                'type':'pic'})
        assert  snp_inst.value == '0.2'
    @staticmethod
    def test_add_contig_annot():
        '''It test the simple case - add contig_annot'''
        annotation_dict = {'name': 'Contig1', 'type':'blast', 'value':'GO:0022'}

        engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_snp_miner_database(engine)
        fhand = open(os.path.join(DATA_DIR, 'example.caf'), 'r')
        add_contig_to_db(engine, fhand)
        add_contig_annot(engine, [annotation_dict])








