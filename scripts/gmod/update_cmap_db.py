#!/usr/bin/env python
''' This scripts updates the feature_acc from cmap's cmap_feature table and
updates it with feature_name '''

from optparse import OptionParser
from franklin.utils.misc_utils import get_db_connection
import os

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()

    parser.add_option('-d', '--db_name', dest='dbname', help='chado db name')
    parser.add_option('-u', '--db_user', dest='dbuser', help='chado db user')
    parser.add_option('-p', '--db_pass', dest='dbpass', help='chado db pass')
    parser.add_option('-t', '--db_host', dest='dbhost', help='chado db host')

    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    dbuser = os.environ['USER'] if options.dbuser is None else options.dbuser
    dbpass = options.dbpass
    dbhost = 'localhost' if options.dbhost is None else options.dbhost
    dbname = 'CMAP' if options.dbname is None else options.dbname

    db_data = {'user' : dbuser, 'name' : dbname, 'pass' : dbpass,
               'host': dbhost, 'adaptor':'mysql'}

    return db_data

def main():
    'The main part'
    # get parameters
    db_data = set_parameters()

    # get a db connection
    connection  = get_db_connection(db_data)

    # get id name mapping
    feat_accs = get_feat_accs(connection)

    #update db changing feature_acc using feature_name in cmap_feature
    update_acc_with_name(connection, feat_accs)

def update_acc_with_name(connection, feat_accs):
    'It updates the acc field with the feature name '
    cursor = connection.cursor()
    for id_, feat_acc in feat_accs:
        statement = "update cmap_feature set feature_acc='%s' where feature_id=%d" %\
                                                                 (feat_acc, id_)
        cursor.execute(statement)

def get_feat_accs(connection):
    'It get s the id name mapping'
    cursor = connection.cursor()
    cursor2 = connection.cursor()
    cursor.execute("select * from cmap_feature")
    feat_accs = {}
    for row in cursor.fetchall():
        cursor2.execute('select map_acc from cmap_map where map_id=%d' % row[2])
        map_name = cursor2.fetchone()[0]
        feat_id   = row[0]
        feat_name = row[4]
        feat_acc = '%s_%s' % (map_name, feat_name)
        if feat_acc not in feat_accs:
            feat_accs[feat_acc] = []
        feat_accs[feat_acc].append(feat_id)

    # now we put A B if there is more than one feature in a mapset
    feat_accs_ = []
    for feat_acc, ids in feat_accs.items():
        if len(ids) > 1:
            for index, id_ in enumerate(ids):
                letter = chr(65 + index)
                final_feat_acc = feat_acc + '_' + letter
                feat_accs_.append((id_, final_feat_acc))
        else:
            feat_accs_.append((ids[0], feat_acc))

    return feat_accs_

if __name__ == '__main__':
    main()


