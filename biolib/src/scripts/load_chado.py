'''
Created on 2009 eka 9

@author: peio
'''
from optparse import OptionParser
from biolib.chado import add_csv_to_chado, connect_database, load_ontology
from biolib.biolib_utils import call
import os

TABLES     = 'db, organism, cv '
ONTOLOGIES = 'cmv_internal_go_lib'



def reload_db_from_scratch(dbname, dbuser, dbpass, dbhost, sql_dump):
    '''It reloads the database from a sql_dump '''
    cmd = ['dropdb', '-U', dbuser,  '-h', dbhost, dbname]
    call(cmd, env={'PGPASS' :dbpass})
    cmd = ['createdb', '-U', dbuser, '-h', dbhost, dbname]
    call(cmd, env={'PGPASS' :dbpass})

    cmd = ['psql', '-U', dbuser, '-h', dbhost, dbname]
    call(cmd, env={'PGPASS' :dbpass}, stdin=sql_dump)


#pylint: disable
def main():
    '''Main section '''
    parser = OptionParser('usage: %prog -d directory [-t tables]', 
                          version='%prog 1.0')
    parser.add_option('-d', '--directory', dest='directory',
                       help='directory where you store the files to load')
    parser.add_option('-t', '--tables', dest='tables',
                       help='Tables to add to chado comma separated')
    parser.add_option('-o', '--ontologies', dest='ontologies',
                       help='own ontologies to addt')
    parser.add_option('-D', '--dbname', dest='dbname', help='chado db name',
                      default= 'chado_test')
    parser.add_option('-u', '--dbuser', dest='dbuser', help='chado db user',
                      default='chado_admin')
    parser.add_option('-p', '--dbpass', dest='dbpass', help='chado db pass',
                      default='chado_pass')
    parser.add_option('-H', '--dbhost', dest='dbhost', help='chado db host',
                      default='localhost')
    
    parser.add_option('-r', '--reload', dest='reload', 
                      help='reload from scratch')
    (options, args) = parser.parse_args()
    

    dbname = options.dbname
    dbuser = options.dbuser
    dbpass = options.dbpass
    dbhost = options.dbhost
    
    if options.directory is None:
        parser.error('Directory with data files is mandatory')
    else:
        directory = options.directory
        
    if options.tables is None:
        tables_ = TABLES
    else:
        tables_ = options.tables
    tables = [table.strip() for table in tables_.split(',')]
    
    if options.ontologies is None:
        ontologies_ = ONTOLOGIES
    else:
        ontologies_ = options.ontologies
    ontologies = [ontology.strip() for ontology in ontologies_.split(',')]
        
        
    if options.reload:
        sql_dump = open(options.reload, 'r')
        reload_db_from_scratch(dbname, dbuser, dbpass, dbhost, sql_dump)
        
    engine = connect_database(dbname, username=dbuser, password=dbpass,
                               host=dbhost)
    
    
    # Let's add info to tables
    no_relation_tables = ['db', 'cv', 'organism']
    for table_name in no_relation_tables:
        if  table_name in tables:
            print "Adding data to %s table" % table_name
            fhand = open('%s/%s.txt' % (directory, table_name) ,'r')
            add_csv_to_chado(fhand, table_name, engine)
    
    # Once I have add cv and  db We can add our ontologies
#    for ontology_file in ontologies:
#        print "Adding ontology %s" % ontology_file
#        fhand_ontology = open(os.path.join(directory,ontology_file), 'rb')
#        load_ontology(fhand_ontology, dbname, dbuser, dbpass, dbhost )
#    
    
     
 
 
 
 
if __name__ == '__main__':
    main()
