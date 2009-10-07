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
from biolib.db.db_utils import  DbMap

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
                                            'rel_attr':'Snv'},
                                  'library_id' :{'kind':'one2many',
                                                 'rel_attr':'library'}}},
                    {'name':'LibrarySnvAlleles',
                     'relations':{'library_snv_id':{'kind':'one2many',
                                                   'rel_attr':'library_svn'}}},
                    {'name':'LibrarySnvProp',
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
            Column('accesion', String, nullable=False, unique=True))
    #reference_table
    Table('Reference', metadata,
          Column('reference_id', Integer, primary_key=True),
          Column('name',   String,  nullable=False, unique=True))

    #reference_annot_table
    Table('ReferenceProp', metadata,
          Column('reference_prop_id', Integer, primary_key=True),
          Column('reference_id', Integer, ForeignKey('Reference.reference_id'),
                 nullable=False),
          Column('type',             String,  nullable=False),
          Column('value',            String,  nullable=False),
          Column('start',            Integer),
          Column('end',              Integer),
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
        Column('kind',   String,  nullable=True),
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
        Column('kind', String, nullable=False),
        Column('num_reads', Integer, nullable=False),
        Column('qualities',  String, nullable=True),
        UniqueConstraint('library_snv_id', 'allele' ))

    #LibrarySnvProp. The annotations associated with an snv in a library
    Table('LibrarySnvProp', metadata,
            Column('library_snv_prop_id', Integer, primary_key=True),
            Column('library_snv_id', Integer,
                   ForeignKey('LibrarySnv.library_snv_id'),  nullable=False),
            Column('kind',  String, nullable=False),
            Column('value', String, nullable=False),
            UniqueConstraint('library_snv_id', 'kind', 'value'))

    metadata.create_all(engine)

def add_snv_to_db(engine, seqvar, library):
    '''It adds the snps to the db'''
    snp_miner = DbMap(engine, SNPMINER_MAP_DEF)
    try:
        # add reference to database
        _insert_reference(snp_miner, seqvar)

        # add location data to database
        _insert_location(snp_miner, seqvar)

        # add library to database
        _insert_library(snp_miner, library)
#        snp_miner.commit()
#        return
        # add seq_var core information
        _insert_seqvar_core(snp_miner, seqvar, library)

        # add alleles to db
        _insert_seqvar_alleles(snp_miner, seqvar, library)

        # add_seqvar_properties to db
        _insert_seqvar_properties(snp_miner, seqvar, library)

        #commit changes
        snp_miner.commit()
    except Exception:
        snp_miner.rollback()
        raise

def _insert_reference(snp_miner, seqvar):
    'It adds the reference data to the database'
    ref_name =  _get_reference_name(seqvar)
    try:
        seq = seqvar.reference.seq
    except AttributeError:
        seq = None

    ref_attr = {'name':ref_name, 'seq':seq}
    snp_miner.get('reference', attributes=ref_attr)

def _insert_location(snp_miner, seqvar):
    'It adds the reference data to the database'
    try:
        name = seqvar.reference.name
    except AttributeError:
        name = seqvar.reference
    position = seqvar.location
    loc_attr = {'reference_id':{'name':name} , 'position':position}
    snp_miner.get('location', attributes=loc_attr)

def _insert_library(snp_miner, library):
    'It adds the library to the database'
    snp_miner.get('library', attributes={'accesion':library})

def _insert_seqvar_core(snp_miner, seqvar, library):
    '''It transform the given parsed file dict in a db dict format'''
    snp_name       = seqvar.name
    ref_name       = seqvar.reference.name
    position       = seqvar.location
    kind           = seqvar.kind

    seqvar_attrs = {'name': snp_name,
                    'location_id': {'reference_id':{'name':ref_name},
                                    'position': position},
                    'type':kind,
                    'library_id':{'accesion':library}}
    snp_miner.get('seqvar', attributes=seqvar_attrs)

def _insert_seqvar_alleles(snp_miner, seqvar, library):
    'Insert alleles of a seq var in database'
    ref_name =  _get_reference_name(seqvar)
    position = seqvar.location
    qual_separator = ":"
    seqvar_name = seqvar.name
    for allele in seqvar.alleles:
        allele_base = allele['allele']
        reads       = allele['reads']
        kind        = allele['kind']
        if 'quality' in allele:
            qual = [str(q) for q in allele['quality']]
            qual    = qual_separator.join(qual)
        else:
            qual    = None
        allele_attrs = {'seqvar_id':{'location_id': {'reference_id':
                                                     {'name':ref_name},
                                                      'position': position},
                                     'library_id':{'accesion':library}},
                        'allele': allele_base, 'type': kind, 'num_reads': reads,
                        'quality':qual}
        snp_miner.get('seqvaralleles', attributes=allele_attrs )

def _insert_seqvar_properties(snp_miner, seqvar, library):
    'Insert alleles of a seq var in database'
    ref_name =  _get_reference_name(seqvar)
    position = seqvar.location
    seqvar_name = seqvar.name

    for annotation, value in seqvar.annotations.items():
        prop_attr = {'seqvar_id':{'location_id': {'reference_id':
                                                              {'name':ref_name},
                                                  'position': position},
                                  'library_id':{'accesion':library}},
                     'type': annotation,
                     'value':value}
        snp_miner.get('seqvarprop', attributes=prop_attr)

def _get_reference_name(seqvar):
    'It return ref base name of a seqvar'
    try:
        name = seqvar.reference.name
    except AttributeError:
        name = seqvar.reference
    return name

def add_reference_annot(engine, list_of_annotation_dicts):
    '''This  function adds annotation information of a contig. Due to the
    variability of the input files, it interfaces accept a list of dict with
    the information to add'''
    for annotation in list_of_annotation_dicts:
        name = annotation['name']
        type_ = annotation['type']
        value = annotation['value']
        if 'start' in annotation:
            start = annotation['start']
        else:
            start = None
        if 'end' in annotation:
            end = annotation['end']
        else:
            end = None
        annot_attr = {'reference_id':{'name':name},
                      'type' : type_,
                      'value': value,
                      'start': start,
                      'end'  : end }
        contig_annot_map = DbMap(engine, SNPMINER_MAP_DEF)
        try:
            contig_annot_map.get('referenceprop', attributes=annot_attr)
            contig_annot_map.commit()
        except Exception:
            contig_annot_map.rollback()
            raise
