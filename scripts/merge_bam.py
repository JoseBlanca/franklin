#!/usr/bin/env python
'''
Thi sscript takes all the bam files from a directory and merges them using
and optional new header.

It adds to each of the alignment on a bam file the reaf group tag, taking it
from the file name. It can be disabled from options

Created on 07/01/2010

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
import os
from tempfile import NamedTemporaryFile
from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.utils.cmd_utils import call
from franklin.sam import (add_header_and_tags_to_sam, merge_sam, sort_bam_sam, bam2sam,
                        sam2bam)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-d', '--work_dir', dest='work_dir',
                    help='input sequence dir')
    parser.add_option('-o', '--outbam', dest='output', help='Output bam')
    parser.add_option('-r', '--reference', dest='reference',
                      help='fasta reference file')
    return parser

def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.work_dir is None:
        parser.error('work_dirwith bams is mandatory')
    else:
        work_dir = options.work_dir

    if options.output is None:
        output = 'out.bam'
    else:
        output = options.output
    if options.reference is None:
        parser.error('Reference is needed')
    else:
        reference = open(options.reference)

    return work_dir, output, reference

def create_header_from_readgroup(readgroups):
    'It creates a bam header from readgroup list'
    header = []
    for readgroup in readgroups:
        head_line = '@RG    ID:%s LB:%s SM:%s' % (readgroup, readgroup,
                                                  readgroup)
        header.append(head_line)
    return  "\n".join(header)

def add_header_and_tags_bams(work_dir, output_dir):
    'it adds readgroupto bams and return added reaadgroups'
    #add to each of the bams the readgroup_tag
    for bam in os.listdir(work_dir):
        if bam.endswith('.bam'):
            #get the readgroup from the name:
            prefix = ".".join(bam.split('.')[:-1])
            sam = open(os.path.join(output_dir, prefix + '.sam'), 'w')

            temp_sam = NamedTemporaryFile(prefix='%s.' % prefix , suffix='.sam')

            bam2sam(os.path.join(work_dir, bam), temp_sam.name)

            add_header_and_tags_to_sam(temp_sam, sam)

            # close and remove temporal stuff
            sam.close()
            temp_sam.close()

def get_opened_sams_from_dir(dir_):
    'It gets all sams from dir'
    sams = []
    for file_ in os.listdir(dir_):
        if file_.endswith('.sam'):
            sams.append(open(os.path.join(dir_, file_)))
    return sams

def main():
    'The script itself'
    #set parameters
    work_dir, output, reference = set_parameters()

    # make a working tempfir
    temp_dir = NamedTemporaryDir()

    # add readgroup tag to each alignment in bam
    add_header_and_tags_bams(work_dir, temp_dir.name)

    # Prepare files to merge
    sams = get_opened_sams_from_dir(temp_dir.name)
    temp_sam = NamedTemporaryFile()

    # merge all the sam in one
    merge_sam(sams, temp_sam, reference)

    # Convert sam into a bam,(Temporary)
    temp_bam = NamedTemporaryFile(suffix='.bam')
    sam2bam(temp_sam.name, temp_bam.name)

    # finally we need to order the bam
    sort_bam_sam(temp_bam.name, output)

    # and make and index of the bam
    call(['samtools', 'index', output], raise_on_error=True)

    temp_dir.close()

if __name__ == '__main__':
    main()
