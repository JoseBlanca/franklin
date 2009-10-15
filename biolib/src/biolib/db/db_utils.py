'''
Created on 2009 eka 18

@author: peio
'''
import sqlalchemy
from sqlalchemy import create_engine, Table
from sqlalchemy.orm import mapper, relation, sessionmaker
from sqlalchemy.orm import exc as orm_exc

def get_db_url(database_conf):
    'It return a usable db url by sqlalchemy'
    if database_conf == 'sqlite':
        return 'sqlite:///%s.db' % database_conf['dbname']
    else:
        return '%s://%s:%s@%s/%s' % (database_conf['dbname'],
                                     database_conf['dbuser'],
                                     database_conf['dbpass'],
                                     database_conf['dbname'] )

def get_foreign_key(table, column):
    '''It returns the table and column names pointed by the foreign key.

    The given column should be a foreign key in the given table.
    The given table and column should be sqlalchemy objects.
    '''
    referer_table    = table
    referer_col      = column
    referenced_table = None
    referenced_col   = None
    for foreign in referer_table.foreign_keys:
        referenced_col = foreign.column
        if referer_col.references(referenced_col):
            referenced_table = foreign.target_fullname.split('.')[0]
            referenced_col = referenced_col.name
            break
    #does the referer table and column really have a foreign key?
    if referenced_col is not None:
        result = referenced_table, referenced_col
    else:
        result = (None, None)
    return result

def _setup_table(table_name, relations, table_classes, row_classes, metadata):
    ''' It setups one table of chado

    It returns a class that represents the Table and a class that represents
    the rows as objects.
    It also requires the already defined table and row classes.
    '''
    #the class that holds the information from the database table
    table = Table(table_name, metadata, autoload=True)
    #we need the names of the columns in the database to check in the
    #row object inits if the given parameters match with the column names
    col_names = [col.name for col in table.columns]
    #if we're defining relations this columns will also be allowable in
    #the inits of the row class
    if relations is not None:
        col_names.extend([rel['rel_attr'] for rel in relations.values()])
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
        for referer_col, relation_ in relations.items():
            relation_name = relation_['rel_attr']
            relation_type = relation_['kind']
            referer_col_obj = table.columns[referer_col]
            referenced_table, referenced_col = get_foreign_key(table,
                                                            referer_col_obj)
            referenced_tbl_obj = table_classes[referenced_table]
            referenced_col_obj = referenced_tbl_obj.columns[referenced_col]
            if relation_type == 'one2many':
                relation_schema[relation_name] = \
                    relation(row_classes[referenced_table],
                       primaryjoin=(referer_col_obj==referenced_col_obj))
            else:
                raise NotImplementedError(\
                   'Only one to many relations are supported at the moment')
        mapper(row, table,  properties=relation_schema)
    table_classes[table_name], row_classes[table_name] = table, row

def setup_mapping(engine, mapping_definitions):
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

    for mapping_definition in mapping_definitions:
        table_name = mapping_definition['name']

        if 'relations' in mapping_definition:
            relations  = mapping_definition['relations']
        else:
            relations = None
        _setup_table(table_name, relations, table_classes, row_classes,
                     metadata)

    return table_classes, row_classes

def connect_database(database, username='', password='', host='localhost'):
    '''It conects to a chado database and returns a engine'''
    db_url ='postgres://%s:%s@%s/%s' % (username, password, host, database)
    return create_engine(db_url)

class DbMap(object):
    'A db orm'

    def __init__(self, engine, mapping_definitions):
        'The init requires an sqlalchemy engine.'
        self._foreign_cache = {} #a ache for the foreign_key function
        self._mapping_definitions = mapping_definitions
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
