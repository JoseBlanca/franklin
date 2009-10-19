'''
Created on 2009 eka 18

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

from sqlalchemy import (Table, Column, Integer, String, MetaData, ForeignKey,
                        UniqueConstraint)
from sqlalchemy.sql import select
from biolib.db.db_utils import  DbMap
from biolib.seqvar.seqvariation import calculate_kind, Snv

SNPMINER_MAP_DEF = [{'name':'Reference'},
                    {'name':'Library'},
                    {'name':'Snv',
                     'relations':{'reference_id':{'kind':'one2many',
                                                  'rel_attr':'reference'}}},
                    {'name':'ReferenceProp',
                     'relations':{'reference_id':{'kind':'one2many',
                                                  'rel_attr':'reference'}}},
                    {'name':'LibrarySnv',
                     'relations':{'snv_id':{'kind':'one2many',
                                            'rel_attr':'snv'},
                                  'library_id' :{'kind':'one2many',
                                                 'rel_attr':'library'}}},
                    {'name':'LibrarySnvAlleles',
                     'relations':{'library_snv_id':{'kind':'one2many',
                                                   'rel_attr':'library_snv'}}},
                    {'name':'LibrarySnvAnnots',
                     'relations':{'library_snv_id':{'kind':'one2many',
                                                    'rel_attr':'library_snv'}}}
                    ]

def create_snp_miner_database(engine):
    '''it creates a snp miner database '''
    metadata = MetaData()
    metadata.bind  = engine
    #pylint: disable-msg=W0612
    #library table
    Table('Library', metadata,
            Column('library_id', Integer, primary_key=True),
            Column('accession', String, nullable=False, unique=True))
    #reference_table
    Table('Reference', metadata,
          Column('reference_id', Integer, primary_key=True),
          Column('name',   String,  nullable=False, unique=True))

    #reference_annot_table
    Table('ReferenceProp', metadata,
          Column('reference_prop_id', Integer, primary_key=True),
          Column('reference_id', Integer, ForeignKey('Reference.reference_id'),
                 nullable=False),
          Column('type', String,  nullable=False),
          Column('value', String,  nullable=False),
          Column('start', Integer),
          Column('end', Integer),
          UniqueConstraint('reference_id', 'type', 'value'))

    #snv table. A location in a reference sequence is an snv
    Table('Snv', metadata,
        Column('snv_id', Integer, primary_key=True),
        Column('reference_id', Integer, ForeignKey('Reference.reference_id'),
               nullable=False),
        Column('location', Integer, nullable=False),
        UniqueConstraint('reference_id', 'location'))

    #LibrarySnv_table. an snv in a particular library
    Table('LibrarySnv', metadata,
        Column('library_snv_id', Integer, primary_key=True),
        Column('name',   String,  nullable=True),
        Column('library_id', Integer, ForeignKey('Library.library_id'),
               nullable=False),
        Column('snv_id',  Integer, ForeignKey('Snv.snv_id'),
               nullable=False),
        UniqueConstraint('library_id', 'snv_id'))

    #LibrarySnvAlleles. The alleles for an snv in a particular library
    Table('LibrarySnvAlleles', metadata,
        Column('library_snv_alleles_id', Integer, primary_key=True),
        Column('library_snv_id', Integer,
               ForeignKey('LibrarySnv.library_snv_id'), nullable=False),
        Column('allele', String, nullable=False),
        Column('kind', Integer, nullable=False),
        Column('reads', Integer, nullable=False),
        Column('qualities',  String, nullable=True),
        UniqueConstraint('library_snv_id', 'allele', 'kind'))

    #LibrarySnvProp. The annotations associated with an snv in a library
    Table('LibrarySnvAnnots', metadata,
            Column('library_snv_annot_id', Integer, primary_key=True),
            Column('library_snv_id', Integer,
                   ForeignKey('LibrarySnv.library_snv_id'),  nullable=False),
            Column('kind',  String, nullable=False),
            Column('value', String, nullable=False),
            UniqueConstraint('library_snv_id', 'kind'))

    metadata.create_all(engine)

class SnvDb(DbMap):
    'It insert, select from the snv database'
    def __init__(self, engine):
        'It initiates the class'
        super(SnvDb, self).__init__(engine, SNPMINER_MAP_DEF)

    def _get_reference(self, reference):
        'It adds the reference data to the database'
        #ref_name =  self._get_reference_name(snv)
        ref_attr = {'name':reference}
        return self.get('Reference', attributes=ref_attr)

    def _get_library(self, library):
        'It adds the library to the database'
        return self.get('Library', attributes={'accession':library})

    def get_snv_sql(self, reference, location):
        'It selects or inserts a row in the svn table'
        reference = self._get_reference(reference)
        loc_attr = {'reference':reference, 'location':location}
        return self.get('Snv', attributes=loc_attr)

    def _get_info_per_library(self, snv_sql, library_info, name=None):
        'It selects or inserts a row in the LibrarySnv table'
        # get library from library_info and map to db
        if 'library' in library_info:
            library_accesion = library_info['library']
            library_sql = self._get_library(library_accesion)
        else:
            library_sql = None
        loc_attr = {'snv':snv_sql, 'library':library_sql, 'name':name}
        return self.get('LibrarySnv', attributes=loc_attr)


    def create_alleles_per_library(self, snv_sql, library_info):
        'It inserts each of the alleles of the info_per_library data object'
        library_snv_sql = self._get_info_per_library(snv_sql, library_info)

        for allele_info in library_info['alleles']:
            allele = allele_info['allele']
            reads  = allele_info['reads']
            kind   = allele_info['kind']
            if 'qualities' in allele_info:
                qual = repr(allele_info['qualities'])
            else:
                qual = None
            attrs = {'library_snv':library_snv_sql, 'allele': allele,
                     'kind': kind, 'reads': reads, 'qualities':qual}
            self.create('LibrarySnvAlleles', attributes=attrs)

    def create_annotations_per_library(self, snv_sql, library_info):
        'It inserts each of the annotations of the info_per_library data object'
        library_snv_sql = self._get_info_per_library(snv_sql, library_info)

        for kind, value in library_info['annotations'].items():
            attrs = {'library_snv':library_snv_sql,
                     'kind':kind,
                     'value':value}
            self.create('LibrarySnvAnnots', attributes=attrs)

    def create_snv(self, snv):
        'it selects or adds a svn to the database'
        reference     = snv.reference
        location      = snv.location
        library_infos = snv.per_lib_info
        snv_sql       = self.get_snv_sql(reference, location)
        for library_info in library_infos:
            self.create_alleles_per_library(snv_sql, library_info)
            self.create_annotations_per_library(snv_sql, library_info)

    def _select_snv_sql(self, reference, location):
        'It selects or inserts a row in the svn table'
        reference = self._get_reference(reference)
        loc_attr = {'reference':reference, 'location':location}
        return self.select_one('Snv', attributes=loc_attr)

    def select_snv(self, reference, location):
        'it returns a snv object giving the reference and the position'
        snv_sql = self._select_snv_sql(reference, location)
        library_snv_sqls = self.select('LibrarySnv', {'snv':snv_sql})
        per_lib_info = []
        for library_snv_sql in library_snv_sqls:
            alleles = []
            for allele in self.select('LibrarySnvAlleles', {'library_snv':
                                                            library_snv_sql}):
                alleles.append({'kind':allele.kind,
                                'reads':allele.reads,
                                'allele':allele.allele,
                                'qualities':eval(allele.qualities)})
            annotations = []

            for annot in self.select('LibrarySnvAnnots', {'library_snv':
                                                          library_snv_sql}):
                annotations.append({'kind':annot.kind,
                                    'value':annot.value})

            per_lib_info.append({'library':library_snv_sql.library.accession,
                                 'alleles':alleles,
                                 'annotations':annotations})
        return Snv(reference=reference, location=location,
                   per_lib_info=per_lib_info)

    def select_snvs(self):
        'An snvs iterator for the database'
        for snv_sql in self.select('Snv'):
            yield self.select_snv(snv_sql.reference.name, snv_sql.location)





