'''
Created on 2009 eka 4

@author: peio
'''
import sqlalchemy
from sqlalchemy import create_engine, Table
from sqlalchemy.orm import mapper, relation
def connect_database(database, username='', password='', host='localhost'):
    '''It conects to a chado database and returns a engine'''
    db_url ='postgres://%s:%s@%s/%s' % (username, password, host, database)
    return create_engine(db_url)
    

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

def add_libraries(session, fhand, library_info):
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
        elif line.startswith('\tname'):
            name = line.split(':')[1]
            library['name'] = name
        elif line.startswith('\ttype'):
            type_ = line.split(':')[1]
            library['type'] = type_
        elif line.startswith('\tcvterms'):
            cvterms = {}
            data_field = line.split(':', 1)[1]
            cvterm_pairs = data_field.split(',')
            for cvterm_pair in cvterm_pairs:
                cvname, cvtermname = cvterm_pair.split(':')
                cvterms[cvname] = cvtermname
            library['cvterms'] = cvterms
        elif line.startswith('\tprops'):
            properties = {}
            data_field = line.split(':', 1)[1]
            property_pairs = data_field.split(',')
            for property_pair in property_pairs:
                cvname, cvtermname = property_pair.split(':')
                properties[cvname] = cvtermname
            library['property'] = properties      
        elif line.startswith('/library_definition'):
            libraries.append(library)
    return libraries  

def add_basic_chado_table(fhand, table):
    '''It adds to chado rows about a given table . It needs a file with this
     format:
         First line: Header. Start with # caracter and after tab delimited 
         colum definition. It have to match with table columns. but not with id
         Next lines: Tab delimited. One line per table row. Values corresponing
          to column definition
    '''
    table_info = _parse_basic_chado_table_file(fhand)
     
def _parse_basic_chado_table_file(fhand):
    '''It does the real parsing, return a dict with  '''
    table_data  = []
    line        = fhand.readline().strip()
    columns     = line[1:].split()
    num_columns = len(columns)
    for column in columns:
        table_data[column] = []
        
    for i, line in enumerate(fhand):
        if line.isspace() or line[0] == '#':
            continue
        items = line.strip().split()
        if len(items) != num_columns:
            raise IndexError('File not well done. Line %d with errors' % i+2)
        for j, item in enumerate(items):
            table_data[columns[j]].append(item)
    return table_data
        
        
     
     
      