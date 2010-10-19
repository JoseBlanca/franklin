#!/usr/bin/env python


# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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
It gets the codon table freq giving a SeqWithQuality file(repr, json, pickle)
'''

from franklin.seq.readers import seqs_in_file
from tempfile import NamedTemporaryFile
from franklin.utils.cmd_utils import call
from optparse import OptionParser

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--inseqfile', dest='inseqfile',
                      help='input sequence file')
    parser.add_option('-o', '--outseqfile', dest='outseqfile',
                      help='output sequence file')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.inseqfile is None:
        parser.error('In file requierd')
    else:
        infhand = open(options.inseqfile)

    if options.outseqfile is None:
        parser.error('out file requierd')
    else:
        outfpath = options.outseqfile

    return infhand, outfpath

def main():
    'the main part'
    # set parameters
    infhand, outfpath = set_parameters()

    # get orfs
    orf_fhand = NamedTemporaryFile()
    for index, orf in enumerate(get_orfs(infhand)):
        orf_fhand.write('>name_%d\n%s\n' % (index, orf))
    orf_fhand.flush()
    print orf_fhand.name
    raw_input()
    # run codon_table_maker
    cmd = ['cusp', '-sequence', orf_fhand.name, '-outfile', outfpath]
    call(cmd)

def get_orfs(infhand):
    'get orfs'
    for seq in seqs_in_file(infhand):
        for orf in seq.get_features(kind='orf'):
            dna = orf.qualifiers['dna']
            pep = orf.qualifiers['pep']

            if check_pep(pep) and check_dna(dna):
                yield dna

def check_pep(pep):
    'it checks that the pep sequence is correct'
    if pep[0] != 'M':
        return False
#    if pep[-1] != '*':
#        return False
    for unknown in ('U', 'O', 'B', 'Z', 'J', 'X'):
        if unknown in pep :
            return False
    return True

def check_dna(dna):
    'It checks that The dna is correct'
    for nucl in dna:
        if nucl.upper() not in ('A', 'T', 'C', 'G'):
            return False
    return True

if __name__ == '__main__':
    main()
