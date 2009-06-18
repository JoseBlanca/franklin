'''
Created on 2009 eka 4

@author: peio
'''
import sqlalchemy
from sqlalchemy.orm import exc as orm_exc
from sqlalchemy.orm import  sessionmaker
from biolib.biolib_utils import call
from biolib.db.db_utils import setup_mapping, get_foreign_key 
from biolib.file_parsers import library_parser
import csv, tempfile   

class Chado(object):
    'A chado orm'
    _mapping_definitions = [
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
    def __init__(self, engine):
        'The init requires an sqlalchemy engine.'
        self._foreign_cache = {} #a ache for the foreign_key function
        self._table_classes, self._row_classes = \
                               setup_mapping(engine, self._mapping_definitions)
        session_klass = sessionmaker(bind=engine)
        self._session = session_klass()

    def create(self, kind, attributes):
        '''It creates an object of the given kind and it adds it to the chado db

        kind - The object kind (e.g. database, feature)
        attributes - a dict with the column names as keys.
        '''
        #which row object?
        row_klass = self._row_classes[kind]

        #now we create the object
        row_inst = row_klass()
        for attr, value in attributes.items():
            setattr(row_inst, attr, value)
        self._session.add(row_inst)
        return row_inst

    def select(self, kind, attributes=None):
        'It returns the query with the given kind with the given attributes'
        if attributes is None:
            attributes = {}
        #do any of the attributes to get an instance from a foreign key?
        attributes = self._resolve_foreign_keys(kind, attributes)
        #which row object?
        row_klass = self._row_classes[kind]

        #pylint: disable-msg=W0142
        #sqlalchemy likes magic
        return self._session.query(row_klass).filter_by(**attributes)

    def select_one(self, kind, attributes):
        'It returns one instance the given kind with the given attributes'
        #which row object?
        return self.select(kind, attributes).one()

    def _get_foreign_key(self, table, column):
        '''It returns the table and column name pointed by the foreign key.
        
        The given column should be a foreign key in the given table.
        '''
        cache_key = table + '%' + column
        if cache_key in self._foreign_cache:
            return self._foreign_cache[cache_key]

        referer_table    = self._table_classes[table]
        referer_col      = referer_table.columns[column]
        result = get_foreign_key(referer_table, referer_col)

        self._foreign_cache[cache_key] = result
        return result

    def _table_index_in_mapping(self, table):
        'It returns the index in the mapping definitions for the given table'
        #we can't use a dict because we need to define the tables in order
        #some tables depend in other
        for index, mapping in enumerate(self._mapping_definitions):
            if table == mapping['name']:
                return index

    def _resolve_foreign_keys(self, table, attributes):
        '''Given a table an a dict with the colnames and values it fixes all
        the foreign keys.
        '''
        new_attrs = {}
        for column, this_col_attributes in attributes.items():
            if isinstance(this_col_attributes, dict):
                #we have to get the instance associated with this dict?
                #we need to know the table and pointed by the foreign key
                #pylint: disable-msg=W0612
                referenced_table, col = self._get_foreign_key(table, column)
                #if we really have a foreign key in this colum (key)
                #we ask for the instance with the value as attributes
                if referenced_table is None:
                    msg  = 'A dict given for a column with no foreign key -> '
                    msg += 'table:%s column:%s dict:%s ' % \
                            (table, column, str(this_col_attributes))
                    raise ValueError(msg)
                try:
                    row_instance = self.select_one(referenced_table,
                                           attributes=this_col_attributes)
                    #from sqlalchemy.orm import exc as orm_exc
                except sqlalchemy.orm.exc.NoResultFound:
                    msg = 'Failed to select a row instance for:\n'
                    msg += '\ttable %s \n\tattributes: %s' % \
                                    (referenced_table, str(this_col_attributes))
                    raise ValueError(msg)
                #where is the table for the current kind defined in the 
                #mapping definitions list?
                indx = self._table_index_in_mapping(table)
                rel_attr = \
                self._mapping_definitions[indx]['relations'][column]['rel_attr']
                new_attrs[rel_attr] = row_instance
            else:
                new_attrs[column] = this_col_attributes
        return new_attrs

    def get(self, kind, attributes):
        '''It returns one attribute with the given attributes and kind
        
        If more than one row match with the unique attributes from the given
        ones an error will be raised.
        If some of the attributes of the non unique ones are changed they will
        be updated.
        '''
        #do any of the attributes to get an instance from a foreign key?
        attributes = self._resolve_foreign_keys(kind, attributes)
        try:
            return self.select_one(kind, attributes)
        except orm_exc.NoResultFound:
            return self.create(kind, attributes)

    def commit(self):
        'It flushes the session and commits all the changes to the database'
        self._session.commit()

    def rollback(self):
        'It rollsback the sesssion'
        self._session.rollback()
################################################################################
# Functions to add chado tables to a chado database   ##########################
################################################################################
def add_csv_to_chado(fhand, table, engine):
    '''It adds to chado rows about a given table.

    It needs a csv file with the field names in the first row, a table name and
    a sqlalchemy engine.
    '''
    chado = Chado(engine)
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
    chado = Chado(engine)
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
    #Naming_schema
    cvname     = parsed_library['cvname']
    cvtermname = parsed_library['cvtermname']
    genus      = parsed_library['genus']
    specie     = parsed_library['specie']
    uniquename = naming.get_uniquename(name=parsed_library['name'],
                                       kind='library')
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
    
    for property_ in parsed_library['properties']:
        items       = property_.split(':')

        cvnamep     = items[0].strip()
        cvtermnamep = items[1].strip()
        value       = items[2].strip()
        uniquename  = naming.get_uniquename(name=parsed_library['name'],
                                       kind='library')
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

