#!/usr/bin/env python

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.
'''
This script is the second step of the sam to snp database pipeline. It takes
some input files and put the information into the database
Created on 01/10/2009

@author: peio
'''
import sqlalchemy
from optparse import OptionParser
from biolib.seqvar.sam_pileup import seqvars_in_sam_pileup
from biolib.seqvar.snp_miner import add_seqvar_to_db
from biolib.db.db_utils import get_db_url

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i sam_pileup, -o req_pos -p pipeline')
    parser.add_option('-i', '--sampileup', dest='samfile',
                      help='Sam pileup file')
    parser.add_option('-r', '--req_pos', dest='req_posfile',
                      help='Required positions file')
    parser.add_option('-f', '--references', dest='references',
                       help='References file')
    parser.add_option('-l', '--library', dest='library',
                       help='Library source name')

    return parser
def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.samfile is None:
        parser.error('Sam file requierd')
    else:
        samfile = open(options.samfile)
    req_posfile = options.req_posfile

    if options.references is None:
        parser.error('Reference file is required')
    else:
        references = open(options.references)
    # Data to insert in the database
    if options.library is None:
        parser.error('Library source name is requires')
    else:
        library = options.library

    #database_conf
    if options.dbhost is None:
        dbhost = 'localhost'
    else:
        dbhost = options.dbhost

    if options.dbname is None:
        parser.error('DB name required')
    else:
        dbname = options.dbname

    if options.dbuser is None:
        parser.error('DN user required')
    else:
        dbuser = options.dbuser

    if options.dbpass is None:
        parser.error('DN password required')
    else:
        dbpass = options.dbpass

    if options.dbengine is None:
        dbengine = 'postgres'
    else:
        dbengine = options.dbengine

    database_conf = {}
    database_conf['dbhost'] = dbhost
    database_conf['dbname'] = dbname
    database_conf['dbuser'] = dbuser
    database_conf['dbpass'] = dbpass
    database_conf['dbengine'] = dbengine

    return samfile, req_posfile, references, library, database_conf

def main():
    'The main part of the script'
    # set parameters
    samfile, req_posfile, references, library, database_conf =\
                                                 set_parameters()
    #get the sevar from requiered base
    seq_vars =  seqvars_in_sam_pileup(samfile, required_positions=req_posfile,
                                      references=references)
    #prepare database conexion:
    db_url = get_db_url(database_conf)
    engine = sqlalchemy.create_engine(db_url)

    #Insert the seqvars to the database
    for seq_var in seq_vars:
        add_seqvar_to_db(engine, seq_var, library=library)

if __name__ == '__main__':
    main()
