#!/usr/bin/env python
'''It calculates a distribution for some sequence property.

It requires a sequence file (optionally with qualities). It calculates the
distribution of some property in the sequence files (like the sequence length).

The distribution types available are:
    - sequence length:        seq_length_distrib
    - quality:                qual_distrib
    - masked sequence length: masked_seq_distrib

The script can generate a plot with the distribution and/or a text file with the
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
import sqlalchemy
from biolib.seqvar.snv_stats import snv_distrib
from biolib.seqvar.seqvariation import snvs_in_file, add_reference_to_svns
from biolib.seqvar.snp_miner import SnvDb

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -f fastafile [-q quality_file]')
    parser.add_option('-i', '--snv_file', dest='snv_file',
                      help='snvs input file')
    parser.add_option('-r', '--ref_file', dest='ref_file',
                      help='reference sequence input fasta')
    parser.add_option('-p', '--plot', dest='plot',
                      help='plot output file')
    parser.add_option('-w', '--distrib', dest='distrib',
                      help='distribution text output file')
    parser.add_option('-t', '--type', dest='kind',
                      help='type of distribution')
    parser.add_option('-d', '--dbname', dest='dbname',
                      help='Database name')
    return parser

#pylint:disable-msg=R0912
def set_parameters():
    '''It set the parameters for this scripts. From de options or from the
     default values'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    io_fhands = {}
    if options.kind is None:
        raise RuntimeError('Distribution type is mandatory (-t distrib_type)')
    else:
        kind = options.kind

    if options.snv_file is None:
        io_fhands['snv_file'] = None
    else:
        io_fhands['snv_file'] = open(options.snv_file, 'r')

    if options.plot is None:
        io_fhands['plot'] = None
    else:
        io_fhands['plot'] = open(options.plot, 'w')

    if options.distrib is None:
        io_fhands['distrib'] = None
    else:
        io_fhands['distrib'] = open(options.distrib, 'w')

    if options.ref_file is None:
        io_fhands['ref_file'] = None
    else:
        io_fhands['ref_file'] = open(options.ref_file, 'r')

    dbname = options.dbname
    return kind, io_fhands, dbname

def main():
    'The main function'
    kind, io_fhands, dbname = set_parameters()
    if dbname is not None:
        engine = sqlalchemy.create_engine('sqlite:///%s' % dbname)
        snvs   = SnvDb(engine).select_snvs()
        if io_fhands['ref_file']:
            snvs   = add_reference_to_svns(snvs, io_fhands['ref_file'])
    else:
        snvs = snvs_in_file(snv_fhand=io_fhands['snv_file'],
                            ref_fhand=io_fhands['ref_file'])
    snv_distrib(snvs, kind=kind, distrib_fhand=io_fhands['distrib'],
                plot_fhand=io_fhands['plot'])

if __name__ == '__main__':
    main()
