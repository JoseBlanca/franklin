'''
Created on 2009 eka 5

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

import unittest
from StringIO import StringIO
import sqlalchemy
from sqlalchemy import (Table, Column, Integer, String, MetaData, ForeignKey,
                        UniqueConstraint)
from sqlalchemy.orm import sessionmaker
from biolib.gmod.chado import (add_csv_to_chado, add_libraries_to_chado,
                               CHADO_MAPPING_DEFINITIONS)
from biolib.db.db_utils import DbMap, setup_mapping
import os, biolib
from tempfile import NamedTemporaryFile
from biolib.db.naming import (create_naming_database,
                                  add_project_to_naming_database,
                                  DbNamingSchema, FileNamingSchema)


DATA_DIR = os.path.join(os.path.split(biolib.__path__[0])[0], 'data')

def create_chado_example():
    '''It creates a chado mini schema with only two tables and one relaion
    in a sqlite memory database '''
    engine = sqlalchemy.create_engine('sqlite:///:memory:')
    metadata = MetaData()
    metadata.bind  = engine
    #aNow I create the tables
    #pylint: disable-msg=W0612
    db_table = Table('db', metadata,
                     Column('db_id', Integer, primary_key=True),
                     Column('name', String, nullable=False, unique=True),
                     Column('description', String))

    dbxref_table = Table('dbxref', metadata,
                         Column('dbxref_id', Integer, primary_key=True),
                         Column('db_id', Integer, ForeignKey('db.db_id'),
                                nullable=False),
                         Column('accession', String, nullable=False),
                         UniqueConstraint('db_id','accession'))
    cv_table = Table('cv', metadata,
                     Column('cv_id', Integer, primary_key=True),
                     Column('name', String, nullable=False, unique=True),
                     Column('definition', String))
    organism_table = Table('organism', metadata,
                     Column('organism_id', Integer, primary_key=True),
                     Column('abbreviation', String),
                     Column('genus', String, nullable=False),
                     Column('species', String, nullable=False),
                     Column('common_name ', String),
                     UniqueConstraint('genus','species'))
    cvterm_table = Table('cvterm', metadata,
                   Column('cvterm_id', Integer, primary_key=True),
                   Column('cv_id', Integer, ForeignKey('cv.cv_id'),
                           nullable=False),
                   Column('name', String, nullable=False),
                   Column('definition', String),
                   Column('dbxref_id', Integer, ForeignKey('dbxref.dbxref_id')),
                          UniqueConstraint('cv_id','name'))
    library_table = Table('library', metadata,
                        Column('library_id', Integer, primary_key=True),
                        Column('organism_id', Integer,
                               ForeignKey('organism.organism_id'),
                               nullable=False),
                        Column('name', String),
                        Column('uniquename', String, nullable=False),
                        Column('type_id', Integer,
                               ForeignKey('cvterm.cvterm_id'), nullable=False),
                       UniqueConstraint('organism_id', 'uniquename', 'type_id'))
    libraryprop_table = Table('libraryprop', metadata,
                            Column('libraryprop_id', Integer, primary_key=True),
                            Column('library_id', Integer,
                                   ForeignKey('library.library_id'),
                                   nullable=False),
                            Column('type_id', Integer,
                               ForeignKey('cvterm.cvterm_id'), nullable=False),
                            Column('value', String),
                            Column('rank', Integer, nullable=False ),
                            UniqueConstraint( 'library_id', 'type_id', 'rank'))
    metadata.create_all(engine)
    return engine

class ChadoOrmTest(unittest.TestCase):
    ''' It test chado orm functions. We use sqlite instead of postgres'''
    @staticmethod
    def test_setup_mapping():
        '''It test basic mapping setup '''
        engine = create_chado_example()
        mapping_definitions = [{'name' :'db'},
                       {'name':'dbxref',
                        'relations':{'db_id':{'kind':'one2many',
                                           'rel_attr':'db'}
                                    }
                       }]
        #pylint: disable-msg=W0612
        table_classes, row_classes = setup_mapping(engine, mapping_definitions)
        #pylint: disable-msg=C0103
        Db       = row_classes['db']
        Dbxref   = row_classes['dbxref']

        Session = sessionmaker(bind=engine)
        session = Session()

        new_db  = Db(name='comav', description='comac_seq')
        session.add(new_db)
        new_dbxref = Dbxref(db=new_db, accession='CMV9789')
        session.add(new_dbxref)
        a_dbxref = session.query(Dbxref).filter_by(accession='CMV9789').first()
        assert isinstance(a_dbxref.dbxref_id, int)
        assert isinstance(a_dbxref.db_id, int)

    def test_chado_orm_helper_basic(self):
        'It test the create, select and get helper functions'
        engine = create_chado_example()
        chado = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
        #######
        #create
        #######
        new_db = chado.create(kind='db', attributes={'name':'hola_db',
                                            'description':'a fake database'})
        assert new_db.name == 'hola_db'
        #if we try to create another one and we flush it we get an error
        try:
            chado.create(kind='db', attributes={'name':'hola_db'})
            chado.commit()
            self.fail()
        except sqlalchemy.exc.IntegrityError:
            chado.rollback()

        #######
        #select
        #######
        new_db = chado.create(kind='db', attributes={'name':'hola_db',
                                            'description':'a fake database'})
        new_dbxref = chado.create(kind='dbxref', attributes={'accession':'666',
                                                             'db':new_db})

        the_db = chado.select_one(kind='db', attributes={'name':'hola_db'})
        assert new_db is the_db
        the_dbxref = chado.select_one(kind='dbxref',
                                      attributes={'accession':'666',
                                                  'db':new_db})
        assert new_dbxref is the_dbxref

        ######
        #get
        ######
        the_db  = chado.get(kind='db', attributes={'name':'hola_db'})
        the_db2 = chado.get(kind='db', attributes={'name':'hola_db'})
        assert the_db is the_db2
        chado.commit()
        chado2 = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
        the_db3  = chado2.select_one(kind='db', attributes={'name':'hola_db'})
        assert the_db3.name == 'hola_db'

    def test_chado_orm_helper_recursive(self):
        '''We can get instances from the session giving dicts inside the id
        columns'''
        engine = create_chado_example()
        chado = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
        #one test with db and dbxref
        chado.get(kind='db', attributes={'name':'hola_db',
                                         'description':'a fake database'})
        chado.get(kind='dbxref', attributes={'accession':'666',
                                             'db_id':{'name':'hola_db'}})
        #a similar example that should fail
        #in this case the get should not create recursively all the instances
        #for all the referenced data.
        try:
            chado.get(kind='dbxref', attributes={'accession':'666',
                                                'db_id':{'name':'hola2_db'}})
            self.fail('ValueError expected')
            #pylint: disable-msg=W0704
        except ValueError:
            pass

EXAMPLE_DBS = \
'''name,description
my_db,some database
another_db,yet another database
'''

class AddCsvChadoTest(unittest.TestCase):
    'It test that we can dump a csv file to a chado database'
    @staticmethod
    def test_add_csv_to_chado():
        'It test that we can dump a csv file to a chado database'
        db_fhand = StringIO(EXAMPLE_DBS)
        engine = create_chado_example()
        #we add the file to the chado table
        add_csv_to_chado(db_fhand, table='db', engine=engine)

        #are they really there?
        chado = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
        db_ins = chado.select_one('db', {'name':'my_db'})
        assert db_ins.description == 'some database'

        #if we try to add it again we won't get more
        dbs = chado.select('db').all()
        assert len(dbs) == 2
        add_csv_to_chado(db_fhand, table='db', engine=engine)
        assert len(dbs) == 2

EXAMPLE_LIBRARY = \
'''format-version:1
library_definition
    name: a
    type: library type:genomic
    organism:Cucumis melo
    cvterms: SO:0001, SO:0002
    properties: property type:strain:Oregon-R, property type:stage:adult male

library_definition
    name:b
    type: library type:genomic
    organism: Cucumis melo
    cvterms:SO:0003, SO:0004
    properties: property type:strain:a_fly, property type:stage:pupa
'''
class AddLibraryToChado(unittest.TestCase):
    'It test that we can dump a library file'
    @staticmethod
    def test_add_libraries():
        'It test that we can add a library file to chado'
        library_fhand = NamedTemporaryFile()
        library_fhand.write(EXAMPLE_LIBRARY)
        library_fhand.flush()
        library_fhand.seek(0)

        engine = create_chado_example()
        chado = DbMap(engine, CHADO_MAPPING_DEFINITIONS)
        # create data necesary to add library
        chado.get(kind='cv', attributes={'name':'library type',
                                         'definition':'a fake cv'})
        chado.get(kind='cv', attributes={'name':'property type',
                                         'definition':'a fake cv'})
        chado.get(kind='db', attributes={'name':'hola_db',
                                         'description':'a fake database'})
        chado.get(kind='organism', attributes={'genus':'Cucumis',
                                               'species':'melo'})
        chado.get(kind='dbxref', attributes={'accession':'001',
                                             'db_id':{'name':'hola_db'}})
        chado.get(kind='dbxref', attributes={'accession':'002',
                                             'db_id':{'name':'hola_db'}})
        chado.get(kind='dbxref', attributes={'accession':'003',
                                             'db_id':{'name':'hola_db'}})
        chado.get(kind='dbxref', attributes={'accession':'004',
                                             'db_id':{'name':'hola_db'}})
        chado.get(kind='cvterm', attributes={'name':'genomic',
                        'cv_id':{'name':'library type'},
                        'dbxref_id':{'accession':'001',
                                     'db_id':{'name':'hola_db'}}})
        chado.get(kind='cvterm', attributes={'name':'strain',
                        'cv_id':{'name':'property type'},
                        'dbxref_id':{'accession':'002',
                                     'db_id':{'name':'hola_db'}}})
        chado.get(kind='cvterm', attributes={'name':'stage',
                        'cv_id':{'name':'property type'},
                        'dbxref_id':{'accession':'003',
                                     'db_id':{'name':'hola_db'}}})

        chado.select_one(kind='cvterm',
                         attributes={'cv_id': {'name': 'library type'},
                                     'name' : 'genomic'})
        #the naming stuff test set up
        fhand = NamedTemporaryFile()
        naming_engine = sqlalchemy.create_engine('sqlite:///:memory:')
        create_naming_database(naming_engine)
        add_project_to_naming_database(naming_engine, name='my_project',
                                       code='my', description='a test project')
        #the real naming stuff
        naming = DbNamingSchema(naming_engine, project='my_project')
        naming = FileNamingSchema(fhand, naming)

        #now we add the libraries to the file
        add_libraries_to_chado(library_fhand, engine, naming)

        library_ins = chado.select_one('library', {'name':'a'})
        assert library_ins
        assert library_ins.type.name == 'genomic'
        assert library_ins.name == 'a'


if __name__ == '__main__':
    unittest.main()
