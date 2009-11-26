#!/usr/bin/env python
'''It calculates the coverage distribution for a samtools pileup file.

It requires a samtools pileup file.The script can generate a plot with the distribution and/or a text file with the
distribution values.
'''

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

from optparse import OptionParser
from biolib.statistics import calculate_read_coverage

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-i', '--pileup', dest='pileup',
                      help='samtools pileup file')
    parser.add_option('-p', '--plot', dest='plot',
                      help='plot output file')
    parser.add_option('-w', '--distrib', dest='distrib',
                      help='distribution text output file')
    parser.add_option('-m', '--max', dest='max', default=None,
                      help='distribution max')
    return parser.parse_args()

#pylint:disable-msg=R0912
def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''

    options = parse_options()[0]

    io_fhands = {}
    if options.pileup is None:
        raise RuntimeError('Need pileup file')
    else:
        io_fhands['pileup'] = open(options.pileup, 'r')

    if options.plot is None:
        io_fhands['plot'] = None
    else:
        io_fhands['plot'] = open(options.plot, 'w')

    if options.distrib is None:
        io_fhands['distrib'] = None
    else:
        io_fhands['distrib'] = open(options.distrib, 'w')
    max_ = options.max
    if max_:
        max_ = int(max_)

    return io_fhands, (None, max_)

def main():
    'The main function'
    io_fhands, range_ = set_parameters()
    if range_ == (None, None):
        range_ = None
    calculate_read_coverage(pileup=io_fhands['pileup'],
                            distrib_fhand=io_fhands['distrib'],
                            plot_fhand=io_fhands['plot'],
                            range_= range_)

if __name__ == '__main__':
    main()
