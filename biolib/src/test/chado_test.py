'''
Created on 2009 eka 5

@author: peio
'''
import unittest
import sqlalchemy
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey
from sqlalchemy.orm import sessionmaker
from biolib.chado_orm import setup_maping


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
 
if __name__ == '__main__':
    unittest.main()