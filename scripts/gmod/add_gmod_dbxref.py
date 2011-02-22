#!/usr/bin/env python

'''
This script adds the dbxref field to the gff. It can add dbxref fields
to point to chado or to point to gbrowse.

To point to chado, data needs to be loaded to chado.

In both cases you need the dbxref database name.

'''

print 'DEPRECATED, TO BE REMOVED ONCE MELONOMICS VERSION 2 IS FINISHED'

import sys, os
from optparse import OptionParser
from franklin.utils.misc_utils import get_db_connection
from franklin.gff import (features_in_gff, get_gff_header, write_gff,
                          add_dbxref_to_feature)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile', help='infile')
    parser.add_option('-o', '--outfile', dest='outfile', help='outfile')
    parser.add_option('-r', '--dbxref', dest='dbxref', help='dbxref db name')
    parser.add_option('-k', '--kind', dest='kind',
                      help='gbrowse2tripal or tripal2gbrowse')
    parser.add_option('-x', '--prefix', dest='prefix',
                      help='prefix used in tripal')

    parser.add_option('-d', '--db_name', dest='dbname', help='chado db name')
    parser.add_option('-u', '--db_user', dest='dbuser', help='chado db user')
    parser.add_option('-p', '--db_pass', dest='dbpass', help='chado db pass')
    parser.add_option('-t', '--db_host', dest='dbhost', help='chado db host')

    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('In file required')
    else:
        infhand = open(options.infile)

    if options.outfile is None:
        outfhand = sys.stdout
    else:
        outfhand =  open(options.outfile, 'w')

    if options.dbxref is None:
        parser.error('Dbxref database required')
    else:
        dbxref = options.dbxref

    kind   = options.kind
    prefix = options.prefix

    dbuser = os.environ['CHADO_DB_USER'] if options.dbuser is None else options.dbuser
    dbpass = os.environ['CHADO_DB_PASS'] if options.dbpass is None else options.dbpass
    dbhost = os.environ['CHADO_DB_HOST'] if options.dbhost is None else options.dbhost
    dbname = os.environ['CHADO_DB_NAME'] if options.dbname is None else options.dbname

    db_data = {'user' : dbuser, 'name' : dbname, 'pass' : dbpass,
               'host': dbhost, 'adaptor':'psql'}

    return infhand, outfhand, dbxref, kind, db_data, prefix

def main():
    'the main part'
    infhand, outfhand, dbxref_db, kind, db_data, prefix = set_parameters()

    # It runs the script
    run(infhand, outfhand, dbxref_db, kind, db_data, prefix)

def run(infhand, outfhand, dbxref_db, kind, db_data, prefix, debug=False):
    'It runs the script'
    # get header
    header = get_gff_header(infhand)

    # get features
    features = features_in_gff(infhand, 3)

    # debugging
    if debug:
        chado_getter = _mock_chado_getter
    else:
        chado_getter = _create_dbxref_id_getter_for_chado(db_data, prefix)

    # add dbxref to features
    dbxref_id_getters = {'gbrowse2tripal': chado_getter,
                         'tripal2gbrowse': _get_dbxref_id_for_gbrowse}

    features = add_dbxref_to_features(features, dbxref_db,
                                      dbxref_id_getters[kind])

    # write features
    write_gff(features, outfhand, header)

def _mock_chado_getter(feature):
    'It simulates a chado getter'
    return 'melo000'


def _create_dbxref_id_getter_for_chado(db_data, prefix):
    def _get_dbxref_id_for_chado(feature):
        '''I t gets the dbxref_id using the prefix and the feature_id that gets from
        chado'''
        db_cursor = get_db_connection(db_data).cursor()
        feature_name = feature['id']
        select = "select * from feature where uniquename='%s'" % feature_name
        #print select
        db_cursor.execute(select)
        row = db_cursor.fetchone()
        if row is None:
            return None
        return prefix + str(row[0])
    return _get_dbxref_id_for_chado

def _get_dbxref_id_for_gbrowse(feature):
    'It gets the dbxref_if to use to link to gbrowse'
    return "%s%%3A%s..%s" % (feature['seqid'],
                             feature['start'], feature['end'])

def add_dbxref_to_features(features, dbxref_db, dbxref_id_getter):
    'It adds dbxref to feature taking into account the kind of dbxref to add'

    for feature in  features:
        dbxref_id  = dbxref_id_getter(feature)
        add_dbxref_to_feature(feature, dbxref_db, dbxref_id)
        yield feature

def test():
    'Test for the script'
    from franklin.utils.misc_utils import DATA_DIR
    from StringIO import StringIO

    # chado to tripal
    gff_in  = open(os.path.join(DATA_DIR, 'map_fis.gff3'))
    gff_out = StringIO()
    dbxref  = 'tripal'
    prefix  = 'melo'
    kind    = 'gbrowse2tripal'
    db_data = None
    run(gff_in, gff_out, dbxref, kind, db_data, prefix, debug=True)
    result = gff_out.getvalue()
    assert 'Dbxref=tripal:melo000;ID=Cm50_A06' in result


    # tripal to chado
    gff_in  = open(os.path.join(DATA_DIR, 'map_fis.gff3'))
    gff_out = StringIO()
    dbxref  = 'melon_gbrowse'
    prefix  = 'melo'
    kind    = 'tripal2gbrowse'
    db_data = None
    run(gff_in, gff_out, dbxref, kind, db_data, prefix, debug=True)
    result = gff_out.getvalue()
    assert 'Dbxref=melon_gbrowse:Chrctg0%3A43286529..43761665' in result

    print "test ok"

if __name__ == '__main__':
    main()
    #test()
