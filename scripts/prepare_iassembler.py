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


'''It prepares a iassembler project using a first mira assembly project.

Using an existing mira assembly project It will create a iassembler project
that can skip the first mira assemply run.
 '''

from optparse import OptionParser
from tempfile import NamedTemporaryFile

import os

IASSEMBLER_INPUT_NAME = 'seqs.fasta'

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-m', '--mira project', dest='mirapath',
                      help='Mira project dir')
    parser.add_option('-i', '--iassembler', dest='iassemblerpath',
                      help='Iassembler project dir')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.mirapath:
        mira_path = os.path.abspath(options.mirapath)
    else:
        parser.error('Mira path mandatory')

    if options.iassemblerpath:
        iassembler_path = os.path.abspath(options.iassemblerpath)
    else:
        parser.error('Mira path mandatory')

    return mira_path, iassembler_path

def get_mira_paths(mira_path):

    mira_dir = os.path.basename(mira_path)
    project_name = mira_dir.split('_')[0]

    unigenes = os.path.join(mira_path, '{0:s}_d_results'.format(project_name),
                            '{0:s}_out.padded.fasta'.format(project_name))
    unigenes_qual = os.path.join(mira_path,
                            '{0:s}_d_results'.format(project_name),
                            '{0:s}_out.padded.fasta.qual'.format(project_name))
    contig_read = os.path.join(mira_path,
                               '{0:s}_d_info'.format(project_name),
                           '{0:s}_info_contigreadlist.txt'.format(project_name))

    return unigenes, unigenes_qual, contig_read

def touch(fname, times = None):
    with file(fname, 'a'):
        os.utime(fname, times)

def process_contig_readlist(fpath_in, fpath_out):
    'It reformats the contig read list of mira in a '
    fhand_in = open(fpath_in)
    fhand_out = open(fpath_out, 'w')
    contig_reads = {}
    for line in fhand_in:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        contig, read = line.split()
        if contig not in contig_reads:
            contig_reads[contig] = []
        contig_reads[contig].append(read)

    for contig, reads in  contig_reads.items():
        contig = 'mira_{0:s}'.format(contig.split('_', 1)[1])
        fhand_out.write(contig + '\t')
        fhand_out.write('\t'.join(reads)+ '\n')
    fhand_out.close()

def main():
    'The main part'

    mira_path, iassembler_path = set_parameters()

    # guess the mira files that we need
    unigenes_fpath, unigenes_qual_fpath, mira_contig_read_fpath = get_mira_paths(mira_path)

    #create the iassembler project dir and subdirs
    if not os.path.exists(iassembler_path):
        os.makedirs(iassembler_path)
    mira_1_dir = os.path.join(iassembler_path,
                             '{0:s}_Assembly'.format(IASSEMBLER_INPUT_NAME),
                             'mira')
    os.makedirs(mira_1_dir)

    # prepare contig readlist for iaasembler
    iassembler_contig_mem_fpath = os.path.join(mira_1_dir, 'CMF10')
    process_contig_readlist(mira_contig_read_fpath, iassembler_contig_mem_fpath)


    # copy files into the iassembler project
    iassembler_unigenes = os.path.join(mira_1_dir, 'mira2.fa')
    iassembler_unigenes_qual = os.path.join(mira_1_dir, 'mira2.fa.qual')
    os.symlink(unigenes_fpath, iassembler_unigenes)
    os.symlink(unigenes_qual_fpath, iassembler_unigenes_qual)

    #create iassembler input files. empty files in iassembler dir because we
    #are going to use the firts mira result as input
    touch(os.path.join(iassembler_path, IASSEMBLER_INPUT_NAME))
    touch(os.path.join(iassembler_path, IASSEMBLER_INPUT_NAME + '.qual'))

    msg  = "To run iassembler you must use this command:\n"
    msg += "iassembler -c -i {0:s}\n".format(IASSEMBLER_INPUT_NAME)
    msg += "From your iassembler dir:{0:s}".format(iassembler_path)
    print msg

if __name__ == '__main__':
    main()