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
This script parses a sam pileup file and get valid sevars (after cleaning and
filtering them).
It writes an evaluable seqvar output file

Created on 30/09/2009

@author: peio
'''
from optparse import OptionParser
import sys
from franklin.snv.sam_pileup import snv_contexts_in_sam_pileup
from franklin.pipelines.pipelines import pipeline_runner


def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i sam_pileup, -o req_pos -p pipeline')
    parser.add_option('-i', '--sampileup', dest='infile',
                      help='list of sam pileup files')
    parser.add_option('-o', '--outfile', dest='outfile',
                      help='Outut file. Required positions file or summary')
    parser.add_option('-p', '--pipeline', dest='pipeline',
                      help='filtering pipeline:snp_basic, snp_exhaustive ')
    parser.add_option('-l', '--libraries', dest='libraries',
                      help='list of libraries')

    return parser

def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Input file requierd')
    else:
        infiles = options.infile.split(',')
        in_fhands = [open(inf) for inf in infiles]

    if options.outfile is None:
        out_fhand = sys.stdout
    else:
        out_fhand = open(options.outfile, 'w')

    if options.libraries is None:
        parser.error('Libraries list is required')
    else:
        libraries = options.libraries.split(',')

    pipeline = options.pipeline
    return in_fhands, out_fhand, pipeline, libraries

def main():
    'The main part of the script'
    in_fhands, out_fhand, pipeline, libraries = set_parameters()

    #get seqvars from sam pileup
    seq_vars_with_context = snv_contexts_in_sam_pileup(in_fhands,
                                                       libraries=libraries)

    #filter/clean seq_vars
    if pipeline is not None:
        seq_vars_with_context = pipeline_runner(pipeline, seq_vars_with_context)
    #remove context to the snv iterator
    #print seq_vars_with_context.next()

    seq_vars = (snv[0] for snv in seq_vars_with_context if snv is not None)

    for snv in seq_vars:
        if snv is not None:
            out_fhand.write(repr(snv))

if __name__ == '__main__':
    main()
