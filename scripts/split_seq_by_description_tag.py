'''
This script takes a seq file and splits the file taking into account one of the
 tags of the description field. It creates as many files as needed
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
import os
from franklin.utils.seqio_utils import seqs_in_file, write_seqs_in_file

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i inputfile -w output_dir')
    parser.add_option('-i', '--infile', dest='infile', help='Input file')
    parser.add_option('-w', '--workdir', dest='workdir', help='Output dir')
    parser.add_option('-t', '--tag', dest='tag', help='tag to use to select')
    parser.add_option('-f', '--fileformat', dest='fileformat',
                       help='file format')
    return parser

def set_parameters():
    '''It sets the parameters for the script.'''
    parser = parse_options()
    options = parser.parse_args()[0]

    if options.infile is None:
        parser.error('Script needs the input file')
    else:
        infhand = open(options.infile)

    if options.workdir is None:
        work_dir = '.'
    else:
        work_dir = options.workdir

    if options.tag is None:
        parser.error('You must specify a tag')
    else:
        tag = options.tag

    if options.fileformat is None:
        format = 'fasta'
    else:
        format = options.fileformat

    return infhand, work_dir, tag, format

def get_item_from_tag(description, tag):
    'get the value of the tag from a string separated by spaces'
    items = description.split(" ")
    for item in items:
        try:
            key, value = item.split(':')
            if key == tag:
                return value
        except ValueError:
            pass

def main():
    'The main function'
    infhand, work_dir, tag, format = set_parameters()
    seqs = seqs_in_file(infhand, format=format)
    tags = {}

    # split seqs by tag. Create a list with all the seqrecords
    for seq in seqs:
        item = get_item_from_tag(seq.description, tag)
        if item not in tags.keys():
            name       = "".join(infhand.name.split('.')[:-1])
            name      += '.' + item + '.' + format
            tags[item] = open(os.path.join(work_dir, name), 'w')

        write_seqs_in_file([seq], tags[item], format=format)
    for files in tags.values():
        files.close()

if __name__ == '__main__':
    main()
