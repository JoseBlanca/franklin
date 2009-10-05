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
This script is the first step of the sam to snp database pipeline. This script
parses a sam pileup file and get valid sevars (after cleaning and filtering
them). It writes a required positions file to use in setp 2

Created on 30/09/2009

@author: peio
'''
from optparse import OptionParser
import sys
from biolib.seqvar.sam_pileup import seqvars_in_sam_pileup, sam_to_required_pos
from biolib.pipelines import pipeline_runner


def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i sam_pileup, -o req_pos -p pipeline')
    parser.add_option('-i', '--sampileup', dest='infile',
                      help='Sam pileup file')
    parser.add_option('-o', '--req_pos', dest='outfile',
                      help='Required positions file')
    parser.add_option('-p', '--pipeline', dest='pipeline', default='snp_clean',
                       help='filtering pipeline')
    return parser

def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Input file requierd')
    else:
        infile = open(options.infile)

    if options.outfile is None:
        outfile = sys.stdout
    else:
        outfile = open(options.outfile, 'w')

    pipeline = options.pipeline

    return infile, outfile, pipeline

def main():
    'The main part of the script'
    sam_pileup, outfile, pipeline = set_parameters()

    #get seqvars from sam pileup
    seq_vars = seqvars_in_sam_pileup(sam_pileup)

    #filter/clean seq_vars
    seq_vars = pipeline_runner(pipeline, seq_vars)

    #From seq_vars to required positions
    sam_to_required_pos(seq_vars, outfile)


if __name__ == '__main__':
    main()