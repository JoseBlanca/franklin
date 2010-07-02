#!/usr/bin/env python
''' This scripts updates the feature_acc from cmap's cmap_feature table and
updates it with feature_name '''

from optparse import OptionParser
from franklin.utils.misc_utils import get_db_connection
import os, getpass
from franklin.gff import features_in_gff

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()

    parser.add_option('-d', '--db_name', dest='dbname', help='chado db name')
    parser.add_option('-u', '--db_user', dest='dbuser', help='chado db user')
    parser.add_option('-p', '--db_pass', dest='dbpass', help='chado db pass')
    parser.add_option('-t', '--db_host', dest='dbhost', help='chado db host')
    parser.add_option('-g', '--gff_file', dest='gff', help='gff file')
    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    dbuser = os.environ['USER'] if options.dbuser is None else options.dbuser
    if options.dbpass is None:
        dbpass = getpass.getpass('Password for user %s: ' % dbuser)
    else:
        dbpass = options.dbpass

    dbhost = 'localhost' if options.dbhost is None else options.dbhost
    dbname = 'CMAP' if options.dbname is None else options.dbname

    db_data = {'user' : dbuser, 'name' : dbname, 'pass' : dbpass,
               'host': dbhost, 'adaptor':'mysql'}
    if options.gff is None:
        parser.error('gff file needed')
    else:
        gff = open(options.gff)

    return db_data, gff

def main():
    'The main part'
    # get parameters
    db_data, gff = set_parameters()

    # get a db connection
    connection  = get_db_connection(db_data)

    # It index the feature accs in the gff by seqid and name
    indexed_feature_accs = index_feat_accs_by_seqid_and_name(gff)

    #update db changing feature_acc using feature_name in cmap_feature
    update_acc_with_name(connection, indexed_feature_accs)

def index_feat_accs_by_seqid_and_name(gff):
    '''It indexes the gff's feature accs  using the seqid and the name'''
    feature_accs = {}
    for feature in features_in_gff(gff, 3):
        seq_id = feature['seqid']
        name   = feature['name']
        index = (seq_id, name)
        if index in  feature_accs:
            continue
        feature_accs[index] = '%s_%s' % index
    return feature_accs

def update_acc_with_name(connection, feature_accs):
    'It updates the acc field with the feature_acc of the gff'
    repeated_features = {}
    cursor = connection.cursor()
    cursor2 = connection.cursor()
    cursor.execute("select * from cmap_feature")
    for row in cursor.fetchall():
        cursor2.execute('select map_acc from cmap_map where map_id=%d' % row[2])
        map_name  = cursor2.fetchone()[0]
        feat_name = row[4]
        feature_acc = feature_accs[(map_name, feat_name)]
        if feature_acc in repeated_features:
            num = repeated_features[feature_acc]
            new_feature_acc = make_unique_feat_accs(feature_acc, num)
            repeated_features[feature_acc] += 1
        else:
            new_feature_acc = feature_acc
            repeated_features[feature_acc] = 65
        statement  = "update cmap_feature set "
        statement += "feature_acc='%s' where feature_id=%d" % (new_feature_acc,
                                                               row[0])

        cursor2.execute(statement)

def make_unique_feat_accs(feature_acc, num):
    'it makes a feature_accs unique giving the number'
    letter = chr(num)
    if letter == 'A':
        feature_acc += '_%s' % letter
    else:
        feature_acc = feature_acc[:-1] + letter
    return feature_acc




if __name__ == '__main__':
    main()


