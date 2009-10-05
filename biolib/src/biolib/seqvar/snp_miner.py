'''
Created on 2009 eka 18

@author: peio
'''
from sqlalchemy import (Table, Column, Integer, String, MetaData, ForeignKey,
                        UniqueConstraint)
from biolib.db.db_utils import  DbMap

SNPMINER_MAP_DEF = [{'name':'reference'},
                    {'name':'library'},
                    {'name':'location',
                     'relations':{'reference_id':{'kind':'one2many',
                                                  'rel_attr':'reference'}}},
                    {'name':'referenceprop',
                     'relations':{'reference_id':{'kind':'one2many',
                                                  'rel_attr':'reference'}}},
                    {'name':'seqvar',
                     'relations':{'location_id':{'kind':'one2many',
                                                 'rel_attr':'location'},
                                  'library_id' :{'kind':'one2many',
                                                 'rel_attr':'library'}}},
                    {'name':'seqvaralleles',
                     'relations':{'seqvar_id':{'kind':'one2many',
                                               'rel_attr':'seqvar'}}},
                    {'name':'seqvarprop',
                     'relations':{'seqvar_id':{'kind':'one2many',
                                           'rel_attr':'seqvar'}}}
                    ]

def create_snp_miner_database(engine):
    '''it creates a snp miner database '''
    metadata = MetaData()
    metadata.bind  = engine
    #pylint: disable-msg=W0612
    reference_table = Table('reference', metadata,
                        Column('reference_id', Integer, primary_key=True),
                        Column('name',   String,  nullable=False, unique=True),
                        Column('seq', String, nullable=True))

    reference_annot_table = Table('referenceprop', metadata,
          Column('reference_prop_id', Integer, primary_key=True),
          Column('reference_id', Integer, ForeignKey('reference.reference_id'),
                 nullable=False),
          Column('type',             String,  nullable=False),
          Column('value',            String,  nullable=False),
          Column('start',            Integer),
          Column('end',              Integer),
          UniqueConstraint('reference_id', 'type', 'value'))

    seqvar_location = Table('location', metadata,
        Column('location_id', Integer, primary_key=True),
        Column('reference_id', Integer, ForeignKey('reference.reference_id'),
               nullable=False),
        Column('position', Integer, nullable=False),
        UniqueConstraint('reference_id', 'position'))

    seqvar_table = Table('seqvar', metadata,
        Column('seqvar_id', Integer, primary_key=True),
        Column('name',   String,  nullable=True),
        Column('type',   String,  nullable=False),
        Column('library_id', Integer, ForeignKey('library.library_id'),
               nullable=False),
        Column('location_id',  Integer, ForeignKey('location.location_id'),
               nullable=False),
        UniqueConstraint('library_id', 'location_id'))

    seqvar_alleles = Table('seqvaralleles', metadata,
        Column('seqvar_alleles_id', Integer, primary_key=True),
        Column('seqvar_id', Integer, ForeignKey('seqvar.seqvar_id'),
               nullable=False),
        Column('allele', String, nullable=False),
        Column('type', String, nullable=False),
        Column('num_reads', Integer, nullable=False),
        Column('quality',  String, nullable=True),
        UniqueConstraint('seqvar_id', 'allele' ))

    seqvar_annot = Table('seqvarprop', metadata,
            Column('seqvarprop_id', Integer, primary_key=True),
            Column('seqvar_id', Integer, ForeignKey('seqvar.seqvar_id'),
                    nullable=False),
            Column('type',  String, nullable=False),
            Column('value', String, nullable=False),
            UniqueConstraint('seqvar_id', 'type', 'value'))

    library = Table('library', metadata,
                    Column('library_id', Integer, primary_key=True),
                    Column('accesion', String, nullable=False, unique=True))

    metadata.create_all(engine)

def add_seqvar_to_db(engine, seqvar, library):
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
        _insert_seqvar_alleles(snp_miner, seqvar)

        # add_seqvar_properties to db
        _insert_seqvar_properties(snp_miner, seqvar)

        #commit changes
        snp_miner.commit()
    except Exception:
        snp_miner.rollback()
        raise

def _insert_reference(snp_miner, seqvar):
    'It adds the reference data to the database'
    try:
        name = seqvar.reference.name
    except AttributeError:
        name = seqvar.reference
    try:
        seq = seqvar.reference.seq
    except AttributeError:
        seq = None

    ref_attr = {'name':name, 'seq':seq}
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

def _insert_seqvar_alleles(snp_miner, seqvar):
    'Insert alleles of a seq var in database'
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
        allele_attrs = {'seqvar_id':{'name': seqvar_name},
                        'allele': allele_base, 'type': kind, 'num_reads': reads,
                        'quality':qual}
        snp_miner.get('seqvaralleles', attributes=allele_attrs )

def _insert_seqvar_properties(snp_miner, seqvar):
    'Insert alleles of a seq var in database'
    seqvar_name = seqvar.name

    for annotation, value in seqvar.annotations.items():
        prop_attr = {'seqvar_id':{'name':seqvar_name}, 'type': annotation,
                     'value':value}
        snp_miner.get('seqvarprop', attributes=prop_attr)


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
