'It gives unique names for the features in the databases'

from sqlalchemy import (Table, Column, Integer, String, Boolean, ForeignKey,
                        MetaData, create_engine, sql)
from db_connection import DbConnection

def _create_naming_database(db_connection):
    'It creates a new empty database to hold the naming schema status'
    #the table definition
    metadata = db_connection.metadata
    Table('last_names', metadata,
         Column('id', Integer, primary_key=True),
         Column('project_id', String(50), ForeignKey("projects.id"),
                nullable=False),
         Column('last_name', String(50), unique=True, nullable=False),
         Column('feature_type', String(50), nullable=False),
         Column('lock', Boolean)
    )
    Table('projects', metadata,
         Column('id', Integer, primary_key=True),
         Column('short_name', String(50), unique=True, nullable=False),
         Column('code', String(2), unique=True, nullable=False),
         Column('description', String(300))
    )
    engine = db_connection.engine
    metadata.create_all(engine)


class _CodeGenerator(object):
    '''This class gives the next code, giving the last one'''
    def __init__(self, seed):
        '''init'''
        self._last_code  = seed
        self._dictionary = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    def next(self):
        '''It return the next code '''
        last_code = self._last_code
        length    = len(last_code)
        letter_part, number_part = '', ''
        next_code = None
        for item  in last_code:
            if item.isdigit():
                number_part += item
            else:
                letter_part += item
        len_letter, len_number = len(letter_part), len(number_part)
        if number_part and number_part == len_number * '9':
            if len_letter == 0:
                next_code = "A".ljust(length, "0")
            else:
                next_code = self._next_letter_code(letter_part).ljust(length,
                                                                          "0")
        elif len_number == 0:
            if letter_part == length * "Z":
                next_code = "0" * (length + 1)
            else:
                next_code = self._next_letter_code(letter_part)
        else:
            letters = str(int(number_part) + 1).rjust(len_number, '0')
            next_code = letter_part + letters
        self._last_code = next_code
        return next_code
    def _next_letter(self, letter):
        '''It returns the next letter one caracter'''
        pos = self._dictionary.index(letter)
        return self._dictionary[pos + 1]
    def _next_letter_code(self, letter):
        '''It return the nextletter code'''
        length = len(letter)
        if letter ==  length * "Z":
            return "A" * (length + 1)
        elif letter[-1] == "Z":
            return self._next_letter_code(letter[:-1]) + "A"
        else:
            return letter[:-1] + self._next_letter(letter[-1])

class NamingSchema(object):
    'This class gives unique names for the objects in the databases'
    def __init__(self, project, feature_kind, db_connection):
        '''It inits the NamingSchema.

        project is the name of the project. Independent names will be created
        for each project.
        feature_kind is the type of object that we want to name,
        like EST or Bacend. It should be one SO term.
        db_connection is a DbConnection instance to the database in which the
        naming schema is stored.
        '''
        self._project = project
        self._kind_so = self._get_so_term(feature_kind)
        self._conn = db_connection
        self._code_gen = None
        self._type_code = None  #two letter feature type code
        self._proj_code = None  #two letter project code
        self._curr_code = None  #the last code given

    @staticmethod
    def _get_so_term(feature_kind):
        'It returns the Sequence Ontology term  for the given feature'
        so = {'EST':'SO:0000345', 'transcribed_cluster':'SO:0001457',
             }
        if not feature_kind in so:
            raise ValueError('Unkown feature_kind, please use so')
        return so[feature_kind]

    def _get_type_code(self):
        'It returns the two letter code for the feature kind.'
        type_code = {'SO:0000345': 'ES', 'transcribed_cluster':'UN'}
        kind = self._kind_so
        if not kind in type_code:
            raise ValueError('No code to the feature kind: ' + kind)
        return type_code[kind]

    def _init_code_gen(self):
        'It looks for the last name yielded and it starts the code generator'
        if self._code_gen is not None:
            #already initialized
            return
        #let's look in the database for the last name
        #this is the first name asked we need the last name stored in the
        #database
        #if there is a name previous name in the database we have a seed
        last_db_name = self._get_last_name_in_db()
        #the kind and project codes from the database for this project and type
        type_code = self._get_type_code()
        proj_code = self._conn.get_id('projects', {'short_name':self._project},
                                      'code')
        #if there is a seed the project and type codes should match
        if last_db_name is not None:
            seed_proj_code = last_db_name[:2]
            if seed_proj_code != proj_code:
                msg = 'No match in project codes in database and script: '
                msg += proj_code + ', ' + seed_proj_code
                raise RuntimeError(msg)
            seed_type_code = last_db_name[2:4]
            if seed_type_code != type_code:
                msg = 'No match in type codes in database and script: '
                msg += type_code + ', ' + seed_type_code
                raise RuntimeError(msg)
        else:
            #this is the first time ever, we start a new one
            last_db_name = proj_code + type_code + '0' * 6
        self._type_code = type_code
        self._proj_code = proj_code
        self._curr_code = last_db_name
        self._code_gen = _CodeGenerator(last_db_name[4:])

    def get_next_name(self):
        'It returns the valid name for the next feature in the project'
        code_gen = self._code_gen
        if code_gen is None:
            #the first time we ask for a name we have to init the last current
            #name
            self._init_code_gen()
            code_gen = self._code_gen

        #now we want the next number
        num_next = code_gen.next()
        next_code = self._proj_code + self._type_code + num_next
        self._curr_code = next_code
        return next_code

    def _get_last_name_in_db(self):
        '''It returns the last name for the given feature'''
        project_id = self._conn.get_id('projects', {'short_name':self._project})
        #is there a previous name for the given feature_kind?
        where = {'project_id': project_id, 'feature_type':self._kind_so}
        try:
            last_name = self._conn.get_id('last_names', where=where,
                                                             column='last_name')
        except ValueError:
            last_name = None
        return last_name

    def insert_project(self, code, description=None):
        'It creates a new project'
        values = {'short_name':self._project, 'code':code}
        if description is not None:
            values['description'] = description
        project_id = self._conn.insert('projects', values)
        if not project_id:
            raise ValueError('There was a problem inserting the project')

    def commit_last_name(self):
        '''It stores the current code in the database.
        
        The next time that someone ask for a code to the database, this one
        will be used as seed.
        '''
        current = self._curr_code  #the last code given
        project_id = self._conn.get_id('projects', {'short_name':self._project})
        feature_type = self._kind_so
        table = self._conn.get_table('last_names')
        try:
            #is there already a last name stored for this project and feature?
            name_id = self._conn.get_id(table, where={'project_id':project_id,
                                    'feature_type':feature_type})
            #so we update it
            update = table.update().where(table.c.id == name_id).values(
                                                            last_name=current)
            update = table.update().values(last_name=current)
            res = self._conn.connection.execute(update)
            res.close()
        except ValueError:
            #there is no name yet, so we insert one
            self._conn.insert(table, {'project_id':project_id,
                                      'feature_type':feature_type,
                                      'last_name':current})




def main():
    conn = DbConnection(database=':memory:', drivername='sqlite')
    _create_naming_database(conn)
    #we wet some codes for a project
    names = NamingSchema('melonomics', 'EST', conn)
    names.insert_project(code='ms', description='The melonomics project')
    assert names.get_next_name() == 'msES000001'
    assert names.get_next_name() == 'msES000002'
    #now we store the last code in the database
    names.commit_last_name()
    names = NamingSchema('melonomics', 'EST', conn)
    assert names.get_next_name() == 'msES000003'
    assert names.get_next_name() == 'msES000004'
    names.commit_last_name()
    names = NamingSchema('melonomics', 'EST', conn)
    assert names.get_next_name() == 'msES000005'
    assert names.get_next_name() == 'msES000006'

if __name__ == '__main__':
    main()
