'''
Created on 2009 eka 4

@author: peio
'''
import sqlalchemy
from sqlalchemy.orm import exc as orm_exc
from biolib.utils.cmd_utils import call
from biolib.gmod.file_parsers import library_parser
from biolib.db.db_utils import DbMap
import csv, tempfile

CHADO_MAPPING_DEFINITIONS = [
                    {'name' :'db'},
                    {'name':'cv'},
                    {'name':'organism'},
                    {'name':'dbxref',
                     'relations' :{'db_id':{'kind':'one2many', 'rel_attr':'db'}}
                    },
                    {'name':'cvterm',
                     'relations' :{'cv_id':{'kind':'one2many',
                                            'rel_attr':'cv'},
                                   'dbxref_id':{'kind':'one2many',
                                             'rel_attr':'dbxref'}
                                  }
                    },
                    {'name':'library',
                    'relations':{'organism_id':{'kind':'one2many',
                                             'rel_attr':'organism'},
                                 'type_id':{'kind':'one2many',
                                           'rel_attr':'type'}
                                }
                   },
                   {'name':'libraryprop',
                    'relations':{'type_id':{'kind':'one2many',
                                           'rel_attr':'type'},
                                 'library_id':{'kind':'one2many',
                                           'rel_attr':'library'}
                                 }
                    }
                   ]

################################################################################
# Functions to add chado tables to a chado database   ##########################
################################################################################
def add_csv_to_chado(fhand, table, engine):
    '''It adds to chado rows about a given table.

    It needs a csv file with the field names in the first row, a table name and
    a sqlalchemy engine.
    '''
    chado = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
    try:
        for row in csv.DictReader(fhand, delimiter=','):
            chado.get(kind=table, attributes=row)

    except (orm_exc.NoResultFound, sqlalchemy.exc.IntegrityError), msg:
        chado.rollback()
        raise RuntimeError('''There have been a error during inserts in table %s
         %s''' % (table, msg))
    chado.commit()
def load_ontology(ontology_fhand, dbname, dbuser, dbpass, dbhost):
    '''It adds a ontology to chado.

    It Uses chado-xml and depends on a gmod instalation. You need gmod installed
     locally. It is very important that the obo file is weel formed
    '''
    fileh = tempfile.NamedTemporaryFile()
    cmd = ['go2chadoxml',  ontology_fhand.name]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('go2chadoxml: ' + stderr)
    else:
        fileh.write(stdout)
    fileh.flush()
    cmd = ['stag-storenode.pl', '-d',
        'dbi:Pg:dbname=%s;host=%s;port=5432' % (dbname, dbhost),
        '--user', dbuser, '--password', dbpass, fileh.name]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('stag-storenode.pl: ' + stderr)
def add_libraries_to_chado(fhand, engine, naming):
    '''This function is used to add libraries to chado using a input file.
    The input file format is as follows:
        library_definition
            name : name
            type : type
            organism: GenusSpecies
            cvterms: cvname:cvtermname, cvname:cvtermname, ...
            props: cvname:cvtermname, cvname:cvtermname, ...
        /library_definition
    The naming instance should be able to provide unique names for the
    libraries
    '''
    chado = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
    for parsed_library in library_parser(fhand):
        try:
            #how do we store the library?
            library_attrs  = _library_to_chado_dict(parsed_library, naming)
            #now we store it
            chado.get(kind = 'library', attributes=library_attrs)
            #now for the library properties
            libraryprops = _libraryprop_to_chado_dict(parsed_library, naming)
            for libraryprops_attr in libraryprops:
                chado.get(kind = 'libraryprop', attributes=libraryprops_attr)
            #and now the cvterms
            chado.commit()
        except Exception:
            chado.rollback()
            raise
def _library_to_chado_dict(parsed_library, naming):
    '''It converts the dictionary that returns the parser in a format that
    our chado_insert_function understand.

    The naming should provide unique names for the libraries.
    '''
    cvname     = parsed_library['cvname']
    cvtermname = parsed_library['cvtermname']
    genus      = parsed_library['genus']
    specie     = parsed_library['specie']
    uniquename = naming.get_uniquename(name=parsed_library['name'])
    library_attrs = {'organism_id':{'genus': genus, 'species':specie},
                     'type_id'    :{'cv_id'    :{'name':cvname},
                                    'name'     : cvtermname},
                     'name'       :parsed_library['name'],
                     'uniquename' :uniquename}

    return library_attrs
def _libraryprop_to_chado_dict(parsed_library, naming):
    '''It returns the chado figure to insert lilbrary properties '''
    cvname     = parsed_library['cvname']
    cvtermname = parsed_library['cvtermname']
    genus      = parsed_library['genus']
    specie     = parsed_library['specie']

    for cvnamep, cvtermnamep, value  in parsed_library['properties']:
        uniquename  = naming.get_uniquename(name=parsed_library['name'])
        libraryprop_attrs = {'library_id':{
                               'organism_id':{'genus': genus, 'species':specie},
                               'uniquename':uniquename,
                               'type_id':{'cv_id':{'name':cvname},
                                          'name' : cvtermname}},
                             'type_id':{'cv_id':{'name':cvnamep},
                                        'name' : cvtermnamep},
                             'value' : value,
                             'rank'  : '0'
                            }
        yield libraryprop_attrs

