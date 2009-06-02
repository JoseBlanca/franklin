'''A helper class for simple database manipulations.'''

import sqlalchemy
from sqlalchemy import sql

class DbConnection(object):
    '''This class collects some sql functions to connect to a db.

    It's a thin convenience wrapper around the sqlalchemy connection.
    '''
    #pylint: disable-msg=R0913
    #I think is more clear with those arguments
    def __init__(self, database, username='', password='',
                 drivername='postgres', host='localhost'):
        'It initialites and it creates a connection to the database'
        if drivername == 'sqlite':
            db_url = '%s:///%s' % (drivername, database)
        else:
            db_url = '%s://%s:%s@%s/%s' % (drivername, username, password, host,
                                       database)
        self._engine = None   #the engine
        self._conn = None     #the connection
        self._trans = None    #a possible transaction going on
        self._metadata = None #the metadata
        self._tables = {}       #it will store the table classes
        self._fk_cache = {}
        self._get_id_cache = {}
        self._connect(db_url)

    def _connect(self, db_url):
        'It returns the connection to the database.'
        self._engine = sqlalchemy.create_engine(db_url)
        self._conn = self._engine.connect()
        #the metadata that stores the db schema
        metadata = sqlalchemy.MetaData()
        metadata.bind  = self._engine
        self._metadata = metadata

    def _get_engine(self):
        'It returns the engine associated to the database'
        return self._engine
    engine = property(_get_engine)

    def _get_connection(self):
        'It returns the connection'
        return self._conn
    connection = property(_get_connection)

    def _get_metadata(self):
        'It returns the metadata'
        return self._metadata
    metadata = property(_get_metadata)

    def transaction_begin(self):
        'It starts a transaction'
        if self._trans is not None:
            self._trans.rollback()
            msg = 'Another transaction was going on'
            raise RuntimeError(msg)
        self._trans = self._conn.begin()

    def transaction_commit(self):
        'It commits the current transaction'
        if self._trans is None:
            msg = "There's no transaction was going on"
            raise RuntimeError(msg)
        self._trans.commit()
        self._trans = None

    def _transaction_rollback(self):
        '''It rolls back the transaction.
        
        It there is no open transaction it won't raise an error.
        '''
        if self._trans is None:
            return
        self._trans.rollback()

    def transaction_rollback(self):
        '''It rolls back the transaction.
        
        If there is no transaction open it will raise an error.
        '''
        if self._trans is None:
            raise RuntimeError("There's no transaction going on")
        self._transaction_rollback(self)

    def close(self):
        '''It closes the ongoing transactions and connections.
        
        If there is a transaction it will be rolledback and an error will be
        raised.
        '''
        if self._trans is not None:
            self._transaction_rollback()
        self._conn.close()

    def __del__(self):
        '''Before saying goodbye we close the transactions and connections'''
        if not self._conn.closed:
            self.close()

    def get_table(self, table):
        '''It loads the metadata that define the table and it returns
        the table class.
        
        It does it using the reflexion mechanism from sqlalchemy.
        '''
        #if the table is not an str is already an sqlalchemy table
        if not isinstance(table, str):
            return table
        #if the table metadata is already table we don't have to do it again
        if table in self._tables:
            return self._tables[table]
        #reflecting tables
        metadata = self.metadata
        try:
            self._tables[table] = sqlalchemy.Table(table, metadata,
                                                   autoload=True)
        except sqlalchemy.exc.OperationalError:
            raise RuntimeError('Database not available')
        except sqlalchemy.exc.NoSuchTableError:
            raise ValueError('No table with the given name')
        return self._tables[table]

    @staticmethod
    def _purge_none_items(a_dict):
        'It deletes the items with a None value from the given dict'
        for key, value in a_dict.items():
            if value is None:
                del a_dict[key]

    @staticmethod
    def _where(table, where):
        '''Given some value for the where it returns an sqlalchemy where
        clause.
        
        The given where can be an int or a dict.
        If it's an int it will be assumed that it's the primary key id.
        If it's a dict it every key and value will be the column and value
        used in the where clause.'''
        columns = table.columns
        where_parts = []
        for key, value in where.items():
            col = columns[key]
            #in sqlalchemy the following is called generative selects
            where_parts.append(col == value)
        if len(where_parts) > 1:
            #pylint: disable-msg=W0142
            #i don't know how to do it without this magic *where_parts
            return sql.and_(*where_parts)
        else:
            return where_parts[0]

    def delete(self, table, where):
        '''It removes a row with the given id from the table.
        
        The where can be an int or a dict. If it's an int it will be used
        as the primary key id. It it's a dict it will be used to build the
        where part of the delete statement.
        It returns the number of deleted rows.
        '''
        table_cls = self.get_table(table)
        where_alche = self._where(table_cls, where)
        stmt = table_cls.delete(where_alche)
        result = self.connection.execute(stmt)
        return result.rowcount

    def get_id(self, table, where, column=None):
        '''It returns the primary key id from a row in the given table.
        
        The where isn a dict with the column names as keys and the
        values as values.
        The given values should select only one row from the table, otherwise
        an error will be raised.
        If a column is given it will be used instead of the primary key column.
        '''
        #we need the db table
        table_cls = self.get_table(table)

        #the columns
        if column is None:
            #the primary key columns
            pk_cols = table_cls.primary_key.columns
            #we only support one primary key for this time
            if not pk_cols or len(pk_cols) > 1:
                raise RuntimeError('Primary key not found or multiple')
            col_id = pk_cols.keys()[0]
        else:
            col_id = table_cls.columns[column]

        #for each column in the values we need a where
        whereclause = self._where(table_cls, where)
        #now we can do the query
        query = sql.select([col_id], whereclause)
        result = self.connection.execute(query)
        #we make sure that there is only one row selected
        res_id = result.fetchone()
        #rowcount wouldn't work with some databases, like sqlite
        if not res_id:
            msg = 'No row for the given values in table ' + table_cls.name
            raise ValueError(msg)
        if result.fetchone() is not None:
            msg = 'More than one row for the given values in table ' + \
                                                            table_cls.name
            raise ValueError(msg)
        result.close()
        res_id = res_id[0]
        return res_id

    def _foreign_key_table(self, table, column):
        '''Given a table and a column name it returns the table referenced
        by the foreign key of that column.
        
        If the column is not a foreign key it returns None'''
        cache_key = table + '%' + column
        if cache_key in self._fk_cache:
            return self._fk_cache[cache_key]
        table_cls = self.get_table(table)
        col = table_cls.columns[column]
        fks = col.foreign_keys
        if not fks or len(fks) > 1:
            return None
        table = fks[0].column.table
        self._fk_cache[cache_key] = table
        return table

    def _look_for_the_values_ids(self, table, values):
        '''It replaces the dicts in the values for the id using the foreign
        keys.
        
        For instance if we have a table user with a field email that is a
        foreign key to the table email and we pass a value like
        {'email':'user@domain.com'} it will replace it with {'email':35} if 35
        is the primary key in the email table.
        '''
        #look for primary keys, if there are dicts as values we should look
        #for the primary keys and substitute them
        for row_dict in values:
            for column, value in row_dict.items():
                #is the value an scalar or a dict, if it's a dict we have to
                #look for the scalar in the table referenced by the foreign
                #key
                if not isinstance(value, dict):
                    continue
                #in which table should we look for these values
                #if this column is a foreign key we'll look there for the
                #values
                ref_table = self._foreign_key_table(table, column)
                #if we get no table with the foreign key we guess that
                #the table should be named after the column
                if not ref_table:
                    ref_table = column
                fk_id = self.get_id(ref_table, value)
                row_dict[column] = fk_id

    def insert(self, table, values):
        '''It inserts the given values in a table row.
        
        The values should be a dict with the column names as keys or a list
        of dicts.
        If there is only one dict to insert the primary key id will be
        returned.
        Some of the values for the dictionaries could be replaced by a dict
        with the data need it to ask for the primary key in other table.
        For instance, we could do.
        values={'user':134} or values={'user':{'name':'Dann'}}
        In the second example we will look for the primary key in the table
        user with the name Dann.
        '''
        if not isinstance(values, list):
            values = [values]

        self._look_for_the_values_ids(table, values)

        table_cls = self.get_table(table)
        if len(values) > 1: #it was a list
            #this executes an executemany from the DBI api
            ins_clause = table_cls.insert()
            result = self.connection.execute(ins_clause, values)
            result.close()
            #in this case the last inserted ids are not supported
            return
        else:
            ins_clause = table_cls.insert(values=values[0])
            result = self.connection.execute(ins_clause)
            result.close()
            return result.last_inserted_ids()[0]

def test_chadosql():
    'It test my chado sql class'
    chado = DbConnection(database='chado', username='chado_admin',
                       password='chado_pass', drivername='postgres',
                       host='fisico')
    print chado.delete('dbxref', {'accession':'ac657'})
    print chado.delete('db', {'name':'hola2'})
    print chado.delete('db', {'name':'hola'})
    print chado.delete('db', {'name':'caracola'})
    print chado.insert(table='db', 
                       values=[{'name':'hola'}, {'name':'caracola'}])
    #print db_id
    #print chado.delete('db', db_id)

    print chado.insert(table='dbxref',
                       values={'db_id':{'name':'hola'}, 'accession':'ac657'})

    chado.get_id('db', {'name':'hola'})
    dbxref_id = chado.get_id('dbxref', {'accession':'ac657'})
    print 'dbxref', dbxref_id

    print 'delete ac657', chado.delete('dbxref', {'accession':'ac657'})
    print chado.delete('db', {'name':'hola2'})
    print chado.delete('db', {'name':'hola'})
    print chado.delete('db', {'name':'caracola'})


def main():
    'It runs everything'
    #import tempfile
    #db_fhand = tempfile.NamedTemporaryFile()
    #db_name = db_fhand.name
    #test_sqlalchemy(db_name)
    #test_sqlsoup(db_name)
    #test_sqlalchemy_sql(db_name)
    #db_fhand.close()

    #test_chadosql()

    
if __name__ == '__main__':
    main()

