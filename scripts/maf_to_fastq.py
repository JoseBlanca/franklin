'''
Created on 26/04/2010

@author: peio
'''

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

from optparse import OptionParser
from Bio.SeqIO.Interfaces import SequenceWriter
import sys

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infile', dest='infile',
                    help='maf file')
    parser.add_option('-o', '--outfile', dest='output', help='Output fastq')
    return parser

def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('infile is mandatory')
    else:
        infhand = open(options.infile)

    if options.output is None:
        outfhand = sys.stdout
    else:
        outfhand = open(options.output, 'w')

    return infhand, outfhand

def remove_gaps(seq, qual):
    'It remove the gaps(*) in the seqs'
    new_seq, new_qual = [], []
    for nseq, nqual in filter(lambda x: x[0] != '*', zip(seq, qual)):
        new_seq.append(nseq)
        new_qual.append(nqual)
    return ''.join(new_seq), ''.join(new_qual)



def main():
    'The main thing'
    #set parameters
    infhand, outfhand =  set_parameters()

    for line in infhand:
        line = line.strip()
        if not line:
            continue
        items = line.split(None, 1)
        if len(items) == 1:
            continue
        key, value = items

        if key == 'CO':
            name = value
        if key == 'CS':
            seq = value
        if key == 'CQ':
            qual = value
            if name is not None and seq is not None:
                seq, qual = remove_gaps(seq, qual)
                if len(seq) == len(qual):
                    outfhand.write('@%s\n%s\n+\n%s\n' % (name, seq, qual))
                else:
                    print "Different length seq and qual, seq: %s" % name
            name, seq, qual = None, None, None
    outfhand.close()
    infhand.close()









if __name__ == '__main__':
    main()
