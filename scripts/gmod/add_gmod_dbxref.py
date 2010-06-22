#!/usr/bin/env python

'''
This script adds the dbxref field to the gff. It can add dbxref fields
to point to chado or to point to gbrowse.

To point to chado, data needs to be loaded to chado.

In both cases you need the dbxref database name.


'''


import sys, os
from optparse import OptionParser
from franklin.utils.misc_utils import get_db_connection
def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile', help='infile')
    parser.add_option('-o', '--outfile', dest='outfile', help='outfile')
    parser.add_option('-r', '--dbxref', dest='dbxref', help='dbxref db name')
    parser.add_option('-k', '--kind', dest='kind', help='kind of dbxref to add')
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
    infhand, outfhand, dbxref, kind, db_data, prefix = set_parameters()

    if kind == 'gbrowse2tripal':
        db_cursor = get_db_connection(db_data).cursor()

    elif kind == 'tripal2gbrowse':
        prefix, db_cursor = None, None

    for line in infhand:
        line = line.strip()
        if not line:
            continue
        if line[0] == '#':
            new_line = line
        else:
            gff3_fields = line.split()
            new_gff3_fields = gff3_fields[:8]

            if len(gff3_fields) > 8:
                gff3_attrs = gff3_fields[8]

                # we need to get th dbxref for this kind of analysis
                dbxref_id  = get_dbxref_id(kind, gff3_fields, db_cursor, prefix)
                if dbxref_id is None:
                    new_gff3_attrs = gff3_attrs
                else:
                    # once we have the dbxref we add it to the gff3 attributes field
                    new_gff3_attrs = add_dbxref(gff3_attrs, dbxref, dbxref_id)

                new_gff3_fields.append(new_gff3_attrs)
            new_line = '\t'.join(new_gff3_fields)
        outfhand.write(new_line + '\n')
    outfhand.close()

def get_dbxref_id(kind, gff3_fields, db_cursor=None, prefix=None):
    'It uses different functions to  guess dbxref_id depending on the kind'
    if kind == 'tripal2gbrowse':
        return get_dbxref_id_for_gbrowse(gff3_fields)
    elif kind == 'gbrowse2tripal':
        return get_dbxref_id_from_chado(gff3_fields[8], prefix, db_cursor)


def get_dbxref_id_for_gbrowse(gff3_fields):
    'It gets the dbxref_if to use to link to gbrowse'
    return "%s%%3A%s..%s" % (gff3_fields[0], gff3_fields[3], gff3_fields[4])

def get_dbxref_id_from_chado(attrs, prefix, db_cursor):
    '''I t gets the dbxref_id using the prefix and the feature_id that gets from
    chado'''
    #first we need a field to look at chado
    unique_name = _get_value_from_gff3_attrs(attrs, field='ID')
    select = "select * from feature where uniquename='%s'" % unique_name
    #print select
    db_cursor.execute(select)
    row = db_cursor.fetchone()
    if row is None:
        return None
    return '%s%d' % (prefix, row[0])


def _get_value_from_gff3_attrs(attrs, field='ID'):
    'It gets the the dbxref from a field of the attrs field. Fallbacks to ID'
    for items in attrs.split(';'):
        key, value = items.split('=')
        if key == field:
            return value
    raise ValueError("the field  %s is not in this gff3 attr field" % field)

def add_dbxref(attrs_, dbxref, dbxref_id):
    'It modifies the attr field adding a dbxref'
    attrs = {}
    for items in attrs_.split(';'):
        key, value = items.split('=')
        attrs[key] = value

    dbxref_field = '%s:%s' % (dbxref, dbxref_id)
    if 'dbxref' in attrs:
        attrs['dbxref'] += ',' + dbxref_field
    else:
        attrs['dbxref'] = dbxref_field

    new_attrs = []
    for key, value in attrs.items():
        new_attrs.append('%s=%s' % (key, value))
    return ';'.join(new_attrs)

if __name__ == '__main__':
    main()
