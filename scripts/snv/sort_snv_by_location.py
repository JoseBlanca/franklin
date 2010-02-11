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
This script sorts the snv output file and get required postions

Created on 14/10/2009

@author: peio
'''
import sys
from optparse import OptionParser

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i snv_file, -o output_file')
    parser.add_option('-i', '--input', dest='infile',
                      help='Snv input filefile')
    parser.add_option('-o', '--output', dest='outfile',
                      help='outputfile')
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
    return infile, outfile

def get_and_sort(fhand):
    'It takes the file and sorts it in a list'
    req_positions = set()
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('reference'):
            items = line.split(',')
            cromosome = items[0].split('=')[1]
            location  = items[1].split('=')[1]
            if cromosome[0] == "'":
                cromosome = cromosome[1:-1]

            req_pos = '%s@@%s' % (cromosome, location)
            req_positions.add(req_pos)
    req_positions = list(req_positions)
    return sorted(req_positions, _req_pos_sort_funct)

def _req_pos_sort_funct(pos1, pos2):
    'sorting function'
    cromosome1, loc1 = pos1.split('@@')
    cromosome2, loc2 = pos2.split('@@')

    if cromosome1 > cromosome2:
        return 1
    elif cromosome1 == cromosome2:
        if int(loc1) >= int(loc2):
            return 1
        elif int(loc1) == int(loc2):
            return 0
        else:
            return -1
    else:
        return -1

def main():
    'The main part'
    infile, outfile = set_parameters()

    #get sorted positions positions
    positions = get_and_sort(infile)
    # write to file
    for req_pos in positions:
        cromosome, locations = req_pos.split('@@')
        outfile.write('%s\t%s\n' % (cromosome, locations))

if __name__ == '__main__':
    main()
