#!/usr/bin/env python
'''
This script changes the names of sequences in a file. It takes the name from a
naming schema database.
If the database does not exists, it creates it.
If the project does not exist it creates it.

Created on 21/12/2009

@author: peio

from optparse import OptionParser
'''
from optparse import OptionParser
import sys, sqlalchemy, os
from biolib.db.naming import (create_naming_database, project_in_database,
                              add_project_to_naming_database, DbNamingSchema,
                              change_names_in_files, FileNamingSchema)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-s', '--infile', dest='infile',
                    help='input sequence file')
    parser.add_option('-o', '--outseqfile', dest='outfile',
                      help='output file')
    parser.add_option('-f', '--fileformat', dest="format",
                      help='input file format', default='fasta')
    parser.add_option('-d', '--database', dest='database',
                    help='path to the naming schema database')
    parser.add_option('-p', '--project', dest='project',
                      help='Project name')
    parser.add_option('-t', '--tag', dest="tag",
                      help='Type of seq to change. EST, CONTIG, ...')
    parser.add_option('-c', '--filecache', dest="filecache",
                  help='If you want to give it a file to add or use as a cache')
    return parser


def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Input file is mandatory')
    else:
        infhand = open(options.infile)

    if options.outfile is None:
        outfhand = sys.stdout
    else:
        outfhand = open(options.outfile, 'w')

    format    = options.format
    database  = options.database
    filecache = options.filecache

    if options.project is None:
        parser.error('Project is mandatory')
    else:
        items = options.project.split(',')
        try:
            description = items[2]
        except IndexError:
            description = None
        try:
            code = items[1]
        except IndexError:
            code = None
        project = {'name':items[0], 'code':code, 'description':description}

    if options.tag is None:
        parser.error('type of tag to change is mandatory')
    else:
        tag = options.tag

    # check if database and filecache are None. One of them is necessary
    if database is None and filecache is None:
        parser.error('Database or filecache is mandatory!')

    return infhand, outfhand, format, database, filecache, project, tag

def main():
    'The main part of the script'
    infhand, outfhand, format, database, filecache, project, tag = \
                                                                set_parameters()

    # check if the database exits
    if database is not None:
        engine   = sqlalchemy.create_engine( 'sqlite:///%s'  % database)
        if not os.path.exists(database):
            create_naming_database(engine)
        # check if the project exists
        project_name = project['name']
        if not project_in_database(engine, project_name):
            if project['code']:
                add_project_to_naming_database(engine, name=project_name,
                                           code=project['code'],
                                           description=project['description'])
            else:
                msg = 'To add a new project to the database the code should be given'
                raise ValueError(msg)

        # create a naming schema
        naming = DbNamingSchema(engine, project_name)
    else:
        naming = None

    if filecache is not None:
        if os.path.exists(filecache):
            mode = 'a'
        else:
            mode = 'w'
        fhand = open(filecache, mode)
        naming = FileNamingSchema(fhand, naming)

    try:
        # change name to seqs in file
        change_names_in_files(infhand, outfhand, naming, format, tag)
        naming.commit()
    except:
        #if we don't commit we lose the changes in the database
        raise

if __name__ == '__main__':
    main()
