#!/usr/bin/env python
'''
Thi sscript takes all the bam files from a directory and merges them using
and optional new header.

It adds to each of the alignment on a bam file the reaf group tag, taking it
from the file name. It can be disabled from options

Created on 07/01/2010

@author: peio
'''

from optparse import OptionParser
import os, re
from biolib.utils.misc_utils import NamedTemporaryDir
from biolib.sam import add_tag_to_bam, merge_bam

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-d', '--work_dir', dest='work_dir',
                    help='input sequence dir')
    parser.add_option('-o', '--outbam', dest='output', help='Output bam')

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

    return work_dir, output

def create_header_from_readgroup(readgroups):
    'It creates a bam header from readgroup list'
    header = []
    for readgroup in readgroups:
        head_line = '@RG    ID:%s LB:%s SM:%s' % (readgroup, readgroup,
                                                  readgroup)
        header.append(head_line)
    return  "\n".join(header)

def add_readgroup_to_bams(work_dir, output_dir):
    'it adds readgroupto bams and return added reaadgroups'
    #add to each of the bams the readgroup_tag
    readgroups = []
    for bam in os.listdir(work_dir):
        if bam.endswith('.bam'):
            #get the readgroup from the name:
            readgroup = re.sub('lib_*', '', bam.split('.')[0])
            readgroups.append(readgroup)
            newbam = os.path.join(output_dir, bam)
            bam    = os.path.join(work_dir, bam)
            add_tag_to_bam(bam, {'RG':readgroup}, newbam)
    return readgroups

def main():
    'The script itself'
    #set parameters
    work_dir, output = set_parameters()

    # make a working tempfir
    temp_dir = NamedTemporaryDir().name

    # add readgroup tag to each alignment in bam and get this readgroups
    readgroups = add_readgroup_to_bams(work_dir, temp_dir)

    # create the header from readgroup catched
    header = create_header_from_readgroup(readgroups)

    # merge all the bam in one
    merge_bam(os.listdir(temp_dir), output, header)


if __name__ == '__main__':
    main()
