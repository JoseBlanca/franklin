#!/usr/bin/env python
'''
It adds redirects to drupal.
It needs to connect to the chado database, get feature_id and feature_name.
Then it construct the redirect using the prefix used in tripal.
'''

from optparse import OptionParser
from franklin.utils.misc_utils import get_db_connection


def parse_options():
    '''Main section '''
    parser = OptionParser('usage: %prog -d directory [-t tables]',
                          version='%prog 1.0')
    parser.add_option('-d', '--chado_dbname', dest='cdbname', help='chado db name',
                      default= 'chado')
    parser.add_option('-D', '--drupal_dbname', dest='ddbname',
                      help='drupal db name', default= 'drupaldb')
    parser.add_option('-u', '--dbuser', dest='dbuser', help='chado db user',
                      default='postgres')
    parser.add_option('-p', '--dbpass', dest='dbpass', help='chado db pass',
                      default=None)
    parser.add_option('-H', '--dbhost', dest='dbhost', help='chado db host',
                      default='localhost')
    parser.add_option('-P', '--prefix', dest='prefix', default='melo',
                      help='prefix used in tripal to construct tripal accesion')

    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    chado_dbname  = options.cdbname
    drupal_dbname = options.ddbname
    dbuser = options.dbuser
    dbpass = options.dbpass
    dbhost = options.dbhost
    prefix = options.prefix

    return chado_dbname, drupal_dbname, dbuser, dbpass, dbhost, prefix

def main():
    'The main part'

    #get parameters
    chado_dbname, drupal_dbname, dbuser, dbpass, dbhost, prefix = \
                                                                set_parameters()

    # get feature_id-feature_name tuples
    db_data = {'user' : dbuser, 'name' : chado_dbname, 'host': dbhost,
               'adaptor':'psql', 'pass':dbpass}

    chado_connection = get_db_connection(db_data)
    feature_id_names = get_feat_id_name(chado_connection)

    # insert the data into drupal_resirec table
    db_data = {'user' : dbuser, 'name' : drupal_dbname, 'host': dbhost,
               'adaptor':'psql', 'pass':dbpass}
    drupal_connection = get_db_connection(db_data)
    drupal_cursor     = drupal_connection.cursor()
    for feature_id, name in feature_id_names:
        if redirect_inserted(name, drupal_cursor):
            continue
        tripal_node = prefix + str(feature_id)
        insert  = "insert into path_redirect(source, redirect, type)"
        insert += " values ('%s', '%s', '301')" % (name, tripal_node)

        drupal_cursor.execute(insert)
    drupal_connection.commit()

def redirect_inserted(name, cursor):
    'Is the redirection already inserted'
    cursor.execute("select * from path_redirect where source='%s'" % name)
    if cursor.rowcount:
        print 'feature %s repeated' % name
    return cursor.rowcount

def get_feat_id_name(connection):
    'It gets the features id and name'
    cursor = connection.cursor()
    cursor.execute("select feature_id, name from feature")
    for id_, name in cursor.fetchall():
        yield (id_, name)

if __name__ == '__main__':
    main()
