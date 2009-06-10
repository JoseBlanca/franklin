'''
Created on 2009 eka 4

@author: peio
'''
import sqlalchemy
from sqlalchemy.orm import exc as orm_exc
from sqlalchemy import create_engine, Table
from sqlalchemy.orm import mapper, relation, sessionmaker
from biolib.biolib_utils import call
import csv, tempfile

def setup_maping(engine, mapping_definitions):
    '''It setups a sqlalchemy mapping for the tables we are interested.

    This function is a way to setup a mapping the the most standard options
    without having to write table classes and row classes for each database
    table.
    It requires a sqlalchemy engine bind to an already setup database and
    the definition of the mapping that we want to set up. 
    The mapping definition is a list of dict. For every table in the database
    that we want to set up there's a dict. This dict should have the 'name'
    of the table in the database and an optional 'relations' list. The relations
    should be tuples with the type of the relation ('one_many' or 'many_one') 
    and the name of the property in the row class created that will hold that
    relation. An example mapping definition could be:  
        [{'name' :'db'}, {'name':'dbxref', 'relations':[('one_many', 'db')]}]
    It means that we want to map the tables db and dbxref and that in the 
    objects that map the rows in dbxref we want a property named db that holds
    a relation with the db objects.
    This function will return two dicts one with the classes that represent the
    tables in the db and other with the classes that represent the rows.
    '''
    
    metadata = sqlalchemy.MetaData()
    metadata.bind  = engine
    
    table_classes = {}
    row_classes   = {}
    
    
    def setup_table(table_name, relations):
        ''' It setups one table of chado
        
        It returns a class that represents the Table and a class that represents
        the rows as objects.
        '''
        #the class that holds the information from the database table
        table = Table(table_name, metadata, autoload=True)
        #we need the names of the columns in the database to check in the
        #row object inits if the given parameters match with the column names 
        col_names = [col.name for col in table.columns]
        #if we're defining relations this columns will also be allowable in
        #the inits of the row class
        if relations is not None:
            col_names.extend([rel[1] for rel in relations])
        #the __init__ for the row class
        def init(self, **kwargs):
            '''It inits an object in the tables '''
            for arg_name, value in kwargs.items():
                if arg_name not in col_names:
                    msg = 'Unknown colum ' + arg_name + ' in table ' +\
                                                                     table.name
                    raise ValueError(msg)
                setattr(self, arg_name, value)
        #the row class
        row   = type(table_name, (object,), {'__init__': init})
        #now we map the Table class and the row class
        if relations is None:
            mapper(row, table)
        else:
            relation_schema = {}
            for relation_ in relations:
                relation_type = relation_[0]
                related_table = relation_[1]
                if relation_type == 'one2many':
                    relation_schema[related_table] = \
                                            relation(row_classes[related_table])
                else:
                    raise NotImplementedError(\
                       'Only one to many relations are supported at the moment')
            mapper(row, table,  properties=relation_schema)
        return table, row

    for mapping_definition in mapping_definitions:
        table_name = mapping_definition['name']
        
        if 'relations' in mapping_definition:
            relations  = mapping_definition['relations']
        else:
            relations = None
        table_classes[table_name], row_classes[table_name] = \
                                            setup_table(table_name, relations) 
    
    return table_classes, row_classes

def connect_database(database, username='', password='', host='localhost'):
    '''It conects to a chado database and returns a engine'''
    db_url ='postgres://%s:%s@%s/%s' % (username, password, host, database)
    return create_engine(db_url)
    
class Chado(object):
    'A chado orm'
    # The mapping for a chado database
    _mapping_definitions = [
                    {'name' :'db'},
                    {'name':'cv'},
                    {'name':'organism'}, 
                    {'name':'dbxref', 'relations' :[('one2many', 'db')]},
                    {'name':'cvterm', 'relations' :[('one2many', 'cv'), 
                                                    ('one2many', 'dbxref')]},
                    {'name':'library', 'relations':[('one2many', 'organism'),
                                                    ('one2many', 'cvterm')]}
                   ]
    def __init__(self, engine):
        'The init requires an sqlalchemy engine.'
        self._table_classes, self._row_classes = \
                               setup_maping(engine, self._mapping_definitions)
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
        #which row object?
        row_klass = self._row_classes[kind]

        #pylint: disable-msg=W0142
        #sqlalchemy likes magic
        return self._session.query(row_klass).filter_by(**attributes)

    def select_one(self, kind, attributes):
        'It returns one instance the given kind with the given attributes'
        #which row object?
        return self.select(kind, attributes).one()

    def get(self, kind, attributes):
        '''It returns one attribute with the given attributes and kind
        
        If more than one row match with the unique attributes from the given
        ones an error will be raised.
        If some of the attributes of the non unique ones are changed they will
        be updated.
        '''
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

def add_libraries(fhand, engine):
    '''This function is used to add libraries to chado using a input file.
    The input file format is as follows:
        library_definition
            name : name
            type : type
            cvterms: cvname:cvtermname, cvname:cvtermname, ...
            props: cvname:cvtermname, cvname:cvtermname, ...
        /library_definition
    '''
    libraries = _parse_libraries_file(fhand)
    chado     = Chado(engine)
    for library in libraries:
#        for cv in cvterms:
#            cv_obj = chado.select_one('cv', name=cv)
#            cv_obj.cv_id
#            for term_name in cv:
#                
#        chado.select_one('cvterm', )
       # print library
       pass
def _parse_libraries_file(fhand):
    ''' It parses the file and returns a list with dictionaries for each 
    library'''
    libraries = []
    for line in fhand:
        if line.isspace() or line[0] == '#':
            continue
        line = line.strip()
        if line.startswith('library_definition'):
            library = {}
        elif line.startswith('name'):
            name = line.split(':')[1]
            library['name'] = name.strip()
        elif line.startswith('type'):
            type_ = line.split(':')[1]
            library['type'] = type_.strip()
        elif line.startswith('cvterms'):
            cvterms = {}
            data_field = line.split(':', 1)[1]
            cvterm_pairs = data_field.split(',')
            for cvterm_pair in cvterm_pairs:
                cvname, cvtermname = cvterm_pair.split(':')
                cvname, cvtermname = cvname.strip(), cvtermname.strip()
                if cvname not in  cvterms:
                    cvterms[cvname] = []
                cvterms[cvname].append(cvtermname)
            library['cvterms'] = cvterms
        elif line.startswith('props'):
            properties = {}
            data_field = line.split(':', 1)[1]
            property_pairs = data_field.split(',')
            for property_pair in property_pairs:
                cvname, cvtermname = property_pair.split(':')
                cvname, cvtermname = cvname.strip(), cvtermname.strip()
                if cvname not in  properties:
                    properties[cvname] = []
                properties[cvname].append(cvtermname)
            library['property'] = properties      
        elif line.startswith('/library_definition'):
            libraries.append(library)
    return libraries  

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

def add_own_ontolgy(ontology_fname, dbname, dbuser, dbpass, dbhost):
    '''It adds a ontology to chado. 
    
    It Uses chado-xml and depends on a gmod instalation. You need gmod installed
     locally. It is very important that the obo file is weel formed'''
    fileh = tempfile.NamedTemporaryFile()
    cmd = ['go2chadoxml',  ontology_fname, '>', fileh.name]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('go2chadoxml: ' + stderr)
    cmd = ['stag-storenode.pl', '-d', 
        'dbi:Pg:dbname=%s;host=%s;port=5432' % (dbname, dbhost),
        '--user', dbuser, '--password', dbpass, fileh.name]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('stag-storenode.pl: ' + stderr)
