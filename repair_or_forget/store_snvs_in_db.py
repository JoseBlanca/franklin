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
This script is the second step of the sam to snp database pipeline. It takes
some input files and put the information into the database
Created on 01/10/2009

@author: peio
'''
import sqlalchemy, os
from optparse import OptionParser
from franklin.snv.sam_pileup import seqvars_in_sam_pileup
from franklin.snv.snp_miner import SnvDb, create_snp_miner_database
from franklin.pipelines import pipeline_runner
from franklin.collections_ import RequiredPosition

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i sam_pileup, -o req_pos -p pipeline')
    parser.add_option('-i', '--sampileup', dest='samfile',
                      help='Sam pileup file')
    parser.add_option('-r', '--req_pos', dest='req_posfile',
                      help='Required positions file')
    parser.add_option('-f', '--references', dest='references',
                       help='References file')
    parser.add_option('-p', '--pipeline', dest='pipeline',
                       help='Another filtering pipeline')
    parser.add_option('-l', '--library', dest='library',
                       help='File library')

    parser.add_option('-d', '--dbname', dest='dbname', help='db name')
#    parser.add_option('-u', '--dbuser', dest='dbuser', help='db user')
#    parser.add_option('-w', '--dbpass', dest='dbpass', help='db pass')
#    parser.add_option('-e', '--dbengine', dest='dbengine', help='db engine')
#    parser.add_option('-h', '--dbhost', dest='dbhost', help='db host')

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
    library     = options.library

    if options.references is None:
        references = None
    else:
        references = open(options.references)
    # Data to insert in the database
    pipeline = options.pipeline

    if options.dbname is None:
        parser.error('DB name required')
    else:
        dbname = options.dbname

    return samfile, req_posfile, references, pipeline, dbname, library

def main():
    'The main part of the script'
    # set parameters
    samfile, req_posfile, references, pipeline, dbname, library = \
                                                                set_parameters()

    req_pos_obj = RequiredPosition(open(req_posfile))
    #get the sevar from requiered base
    snvs_with_context =  seqvars_in_sam_pileup(samfile,
                                               required_positions=req_pos_obj,
                                               references=references,
                                               library=library)
    # filter using pipeline
    if pipeline is not None:
        snvs_with_context = pipeline_runner(pipeline, snvs_with_context)

    # remove context
    snvs = (snv[0] for snv in snvs_with_context if snv is not None)

    #prepare database conexion:
    db_url = 'sqlite:///%s'  % dbname
    engine = sqlalchemy.create_engine(db_url)
    if not os.path.exists(dbname):
        create_snp_miner_database(engine)
    snv_miner = SnvDb(engine)

    #Insert the seqvars to the database
    for snv in snvs:
        snv_miner.create_snv(snv)
        snv_miner.commit()

if __name__ == '__main__':
    main()
