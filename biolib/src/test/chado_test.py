'''
Created on 2009 eka 5

@author: peio
'''
import unittest
from StringIO import StringIO
import sqlalchemy
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey
from sqlalchemy.orm import sessionmaker
from biolib.chado import setup_maping, Chado, add_csv_to_chado


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
                                nullable=False, unique=True),
                         Column('accession', String, nullable=False, 
                                unique=True))
    metadata.create_all(engine)
    return engine

class ChadoOrmTest(unittest.TestCase):
    ''' It test chado orm functions. We use sqlite instead of postgres'''
    @staticmethod
    def test_setup_mapping():
        '''It test basic mapping setup '''
        engine = create_chado_example()
        mapping_definitions = [{'name' :'db'}, 
                       {'name':'dbxref', 'relations':[('one2many', 'db')]}]
        #pylint: disable-msg=W0612
        table_classes, row_classes = setup_maping(engine, mapping_definitions)
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
        assert (1, 1) == (a_dbxref.dbxref_id, a_dbxref.db_id) 
 
    def test_chado_orm_helper(self):
        'It test the create, select and get helper functions'
        engine = create_chado_example()
        chado = Chado(engine)
        #######
        #create
        #######
        new_db = chado.create(kind='db', attributes={'name':'hola_db',
                                            'description':'a fake database'})
        assert new_db.name == 'hola_db'
        #if we try to create another one and we flush it we get an error
        try:
            chado.create(kind='db', attributes={'name':'hola_db'})
            chado.flush()
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
        the_db = chado.get(kind='db', attributes={'name':'hola_db'})


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
        chado = Chado(engine)
        db_ins = chado.select_one('db', {'name':'my_db'})
        assert db_ins.description == 'some database'

        #if we try to add it again we won't get more
        dbs = chado.select('db').all()
        assert len(dbs) == 2
        add_csv_to_chado(db_fhand, table='db', engine=engine)
        assert len(dbs) == 2




if __name__ == '__main__':
    unittest.main()
