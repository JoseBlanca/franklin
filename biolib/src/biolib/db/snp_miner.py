'''
Created on 2009 eka 18

@author: peio
'''
from sqlalchemy import (Table, Column, Integer, String, MetaData, ForeignKey, 
                        UniqueConstraint)
from biolib.db.db_utils import  DbMap
from biolib.file_parsers import snp_summary_parser
from biolib.contig_io import  get_parser
from biolib.read_source import get_read_strain

SNPMINER_MAP_DEF = [{'name':'contig'},
                    {'name':'contigprop',
                     'relations':{'contig_id':{'kind':'one2many',
                                            'rel_attr':'contig'}}},
                    {'name':'snp',
                     'relations':{'contig_id':{'kind':'one2many',
                                               'rel_attr':'contig'}}},
                    {'name':'snpalleles',
                     'relations':{'snp_id':{'kind':'one2many',
                                           'rel_attr':'snp'}}},
                    {'name':'snpprop',
                     'relations':{'snp_id':{'kind':'one2many',
                                           'rel_attr':'snp'}}}
                    ]

def create_snp_miner_database(engine):
    '''it creates a snp miner database ''' 
    metadata = MetaData()
    metadata.bind  = engine
    #pylint: disable-msg=W0612
    contig_table = Table('contig', metadata,
                Column('contig_id', Integer, primary_key=True),
                Column('name',   String,  nullable=False))
    contig_annot_table = Table('contigprop', metadata,
                Column('contig_prop_id', Integer, primary_key=True),
                Column('contig_id',   String, ForeignKey('contig.contig_id'),
                         nullable=False),
                Column('type',             String,  nullable=False),
                Column('value',            String,  nullable=False),
                Column('start',            Integer),
                Column('end',              Integer),
                UniqueConstraint('contig_id', 'type', 'value'))
    snp_table = Table('snp', metadata,
        Column('snp_id', Integer, primary_key=True),
        Column('name',   String,  nullable=False, unique=True),
        Column('type',   String,  nullable=False),
        Column('contig_id', String, ForeignKey('contig.contig_id'),
                nullable=False),
        
        Column('start',  Integer, nullable=False),
        Column('end',    Integer, nullable=False))
    snp_alleles = Table('snpalleles', metadata,
        Column('snpprop_id', Integer, primary_key=True),
        Column('snp_id', Integer, ForeignKey('snp.snp_id'), nullable=False),
        Column('read', String, nullable=False),
        Column('allele', String, nullable=False),
        Column('library',  String, nullable=False),
        Column('strain',  String, nullable=False),
        UniqueConstraint('snp_id', 'read', 'allele' ) )
    snp_annot = Table('snpprop', metadata,
            Column('snpannot_id', Integer, primary_key=True),
            Column('snp_id', Integer, ForeignKey('snp.snp_id'), nullable=False),
            Column('type',  String, nullable=False),
            Column('value', String, nullable=False),
            UniqueConstraint('snp_id', 'type', 'value')
            )   
    metadata.create_all(engine)
def add_contig_to_db(engine, fhand):
    '''It adds contigs to the database giving a caf or ace file '''
    format_    = fhand.name[-3:]
    parser     = get_parser(fhand, format_)
    contig_map = DbMap(engine, SNPMINER_MAP_DEF)
    for name in parser.contig_names():
        try:
            contig_attributes = {'name':name}
            contig_map.get('contig', attributes=contig_attributes) 
            contig_map.commit()
        
        except Exception:
            contig_map.rollback()
            raise
def add_contig_annot(engine, list_of_annotation_dicts):
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
        annot_attr = {'contig_id':{'name':name},
                      'type' : type_,
                      'value': value,
                      'start': start,
                      'end'  : end }
        contig_annot_map = DbMap(engine, SNPMINER_MAP_DEF)
        try:
            contig_annot_map.get('contigprop', attributes=annot_attr)
            contig_annot_map.commit()
        except Exception:
            contig_annot_map.rollback()
            raise
    
          
def add_snp_to_db(engine, fhand, read_source, fhands_libraries):
    '''It adds the snps to the db'''
    snp_miner = DbMap(engine, SNPMINER_MAP_DEF)
    for snp_dict in snp_summary_parser(fhand):
        try:
            snp_attrs = _snp_summary_to_db_dict(snp_dict)
            snp_miner.get('snp', attributes=snp_attrs)
            # It adds props to the snp
            for snp_annot_attrs in _snp_summary_prop_to_db_dict(snp_dict):
                snp_miner.get('snpprop', attributes=snp_annot_attrs)
            # It adds alleles/read info to the database
            for snp_alleles_attr in _snp_summary_alleles_to_db_dict(snp_dict, 
                                                            read_source,
                                                            fhands_libraries):
                snp_miner.get('snpalleles', attributes=snp_alleles_attr)
            snp_miner.commit()
        except Exception:
            snp_miner.rollback()
            raise     




def  _snp_summary_to_db_dict(snp_dict):
    '''It transfrom the given parsed file dict in a db dict format'''
    snp_name    = snp_dict['name'] 
    contig_name = snp_dict['contig']
    start       = snp_dict['start']
    end         = snp_dict['end']
    kind        = snp_dict['kind']
    
    snp_attrs = {'name': snp_name,
                 'contig_id': {'name':contig_name},
                 'start' : start,
                 'end': end,
                 'type':kind
                 }
    return snp_attrs

def _snp_summary_prop_to_db_dict(snp_dict):
    '''It converts to db dict the annotations of the snps'''
    
    snp_name = snp_dict['name']   
    for type_, value in snp_dict['annotations'].items():
        snp_annot_attrs = {'snp_id':{'name': snp_name},
                           'type'  : type_,
                           'value' : value }
        yield snp_annot_attrs

def _snp_summary_alleles_to_db_dict(snp_dict, read_source, fhands_libraries):
    '''It converts to db dict the allele distributions of the snps '''
    snp_name    = snp_dict['name']
    for allele, reads in snp_dict['alleles'].items():
        for read in reads:
            library = read_source.get_library(read)
            strain  = get_read_strain(library, fhands_libraries)
            snp_allele_attrs = {'snp_id'     : {'name': snp_name},
                                'read'       : read, 
                                'allele'     : allele ,
                                'library'    : library ,
                                'strain'     : strain}
            yield snp_allele_attrs
