#!/usr/bin/env python

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.
'''
This script configures a existing empy database to be used as our snp database

Created on 01/10/2009

@author: peio
'''
from optparse import OptionParser
import sqlalchemy
from franklin.seqvar.snp_miner import create_snp_miner_database
from franklin.db.db_utils import get_db_url

def parse_options():
    'It parses the command line arguments'
    optmsg  = 'usage: %prog -d dbname -u user -p password -e dbengine'
    parser = OptionParser(optmsg)
    parser.add_option('-d', '--dbname', dest='dbname', help='db name')
    parser.add_option('-u', '--dbuser', dest='dbuser', help='db user')
    parser.add_option('-p', '--dbpass', dest='dbpass', help='db pass')
    parser.add_option('-p', '--dbengine', dest='dbengine', help='db engine')
    parser.add_option('-p', '--dbhost', dest='dbhost', help='db host')
    return parser

def set_parameters():
    '''It sets the parameters for the script.'''
    parser  = parse_options()
    options = parser.parse_args()[0]

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
    return database_conf


def main():
    'The main part of the script'
    # set parameters
    database_conf = set_parameters()

    # get engine
    db_url = get_db_url(database_conf)
    engine = sqlalchemy.create_engine(db_url)

    #configure seq var database_conf
    create_snp_miner_database(engine)

if __name__ == '__main__':
    main()



