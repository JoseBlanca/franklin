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
    msg = 'Which action are you going to do? change_name, delete, show_names '
    parser.add_option('-A', '--action', dest='action', default='change_names',
                      help=msg)

    parser.add_option('-s', '--infile', dest='infile',
                    help='input sequence file')
    parser.add_option('-o', '--outseqfile', dest='outfile',
                      help='output file')
    parser.add_option('-f', '--fileformat', dest="format",
                      help='input file format', default='fasta')
    parser.add_option('-d', '--database', dest='database',
                    help='path to the naming schema database')
    parser.add_option('-p', '--projectname', dest='projectname',
                      help='Project name')
    parser.add_option('-a', '--actiondesc', dest='actiondesc',
                      help='Action description')
    parser.add_option('-c', '--projectcode', dest='projectcode',
                      help='Project name')
    parser.add_option('-e', '--projectdescription', dest='projectdesc',
                      help='Project description')
    parser.add_option('-k', '--feature_kind', dest="feature_kind",
                      help='Type of seq to change. EST, CONTIG, ...')
    parser.add_option('-l', '--filecache', dest="filecache",
                  help='If you want to give it a file to add or use as a cache')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    action = options.action
    if action not in ('change_name', 'show_names', 'delete'):
        parser.error('Wrong action option')

    io_options   = _parse_io_options(options, parser, action)
    project      = _parse_project_option(options, parser, action)
    feature_kind = _parse_feature_kind_option(options, parser, action)
    database     = _parse_database_option(options, parser, action)
    filecache    = options.filecache
    description  = options.actiondesc

    # check if database and filecache are None. One of them is necessary
    if database is None and filecache is None:
        parser.error('Database or filecache is mandatory!')

    return (io_options, database, filecache, project, feature_kind, action,
            description)

def _parse_io_options(options, parser, action):
    'parse IO options'
    if action == 'change_name':
        if options.infile is None:
            parser.error('Input file is mandatory')
        else:
            infhand = open(options.infile)

        if options.outfile is None:
            outfhand = sys.stdout
        else:
            outfhand = open(options.outfile, 'w')

        format    = options.format
    else:
        infhand, outfhand, format = None, None, None
    return {'infhand':infhand, 'outfhand':outfhand, 'format':format}

def _parse_project_option(options, parser, action):
    'It parses project option '
    # Project_name
    if action in ('change_names', 'delete') and options.projectname is None:
        parser.error('Project is mandatory')
    else:
        project_name = options.projectname
    if action == 'change_name':
        if options.projectcode is None:
            parser.error('project Code is mndatory for action: %s' % action)
        else:
            code = options.projectcode
        description = options.projectdesc
    else:
        code, description = None, None
    return  {'name':project_name, 'code':code, 'description':description}

def _parse_feature_kind_option(options, parser, action):
    'It parses feature_kind option'
    if action  not in ('show_names'):
        if options.feature_kind is None:
            parser.error('Type of feature_kind to change is mandatory')
        else:
            feature_kind = options.feature_kind
    else:
        feature_kind = None
    return feature_kind

def _parse_database_option(options, parser, action):
    'It parses database option'
    if action == 'change_name':
        database = options.database
    elif options.database is None:
        parser.error('type of tag to change is mandatory')
    else:
        database = options.database
    return database

def main():
    'The main part of the script'
    (io_options, database, filecache, project, feature_kind,
                                        action, description) = set_parameters()

    if action == 'change_name':
        _change_names_in_file(io_options, database, filecache, project,
                              feature_kind, description)

    elif action == 'delete':
        _delete_last_row(database, project['name'], feature_kind)
        print "The resulting databases is as it is:"
        print _get_names_from_db(database, project_name=project['name'])

    elif action == 'show_names':
        print _get_names_from_db(database, project['name'])

def _get_names_from_db(database, project_name=None):
    'It show the names in the database'
    engine   = sqlalchemy.create_engine( 'sqlite:///%s'  % database)
    naming = DbNamingSchema(engine, project_name, feature_kind=None)
    toprint  = " Project | Project Code |    Name    | Feature | date              | description \n"
    toprint += "-------------------------------------------------------------\n"
    for name in naming.get_names_from_db():
        date = name['date']
        toprint += ' %s |     %s     | %s |   %s  |%d/%d/%d %d:%d:%d | %s\n' % \
                (name['project'], name['project_code'], name['name'],
                 name['feature_type'], date.year, date.month, date.day,
                 date.hour, date.minute, date.second, name['description'])
    return toprint

def _delete_last_row(database, project_name, feature_kind):
    'It deletes the'
    engine   = sqlalchemy.create_engine( 'sqlite:///%s'  % database)
    naming = DbNamingSchema(engine, project=project_name,
                            feature_kind=feature_kind)
    naming.revert_last_name()


def _change_names_in_file(io_options, database, filecache, project,
                          feature_kind, description):
    'It changes the name to a file giving a naming'
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
                msg  = 'To add a new project to the database the code should be'
                msg += ' given'
                raise ValueError(msg)

        # create a naming schema
        naming = DbNamingSchema(engine, project=project_name,
                                feature_kind=feature_kind)
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
        infhand  = io_options['infhand']
        outfhand = io_options['outfhand']
        format   = io_options['format']
        change_names_in_files(infhand, outfhand, naming, format)
        naming.commit(description=description)
    except:
        #if we don't commit we lose the changes in the database
        raise

if __name__ == '__main__':
    main()
