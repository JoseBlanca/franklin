#!/usr/bin/env python

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
'''
It calculates de sam pileup file using one or various bam files.

Created on 08/10/2009

@author: peio
'''

from biolib.biolib_cmd_utils import NamedTemporaryDir, call
import os, sys
from optparse import OptionParser

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i bam_file[,bam_file2] -o sam_pileup')
    parser.add_option('-i', '--bamfiles', dest='bamfiles',
                      help='bam files coma separated, without space')
    parser.add_option('-o', '--sam_pileup', dest='pileup',
                      help='soutput sam_pileup')
    parser.add_option('-r', '--references', dest='reffile',
                      help='reference seq file')
    return parser
def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.bamfiles is None:
        parser.error('bam input file is needed')
    else:
        bamfiles = options.bamfiles.split(',')

    if options.reffile is None:
        parser.error('file is mandatory')
    else:
        reffile = options.reffile

    if options.pileup is None:
        pileup = sys.stdout
    else:
        pileup = open(options.pileup, 'w')

    return bamfiles, pileup, reffile

def main():
    "The real script"
    #set parameters
    bamfiles, pileup, reffile = set_parameters()

    # TEMPORARY DIR NAME
    temp_dir = NamedTemporaryDir()
    dir_name = temp_dir.name

    #merge bam files
    if len(bamfiles) == 1:
        merged_bam_file = bamfiles[0]
    else:
        merged_bam_file = os.path.join(dir_name, 'merged.bam')
        cmd = ['samtools', 'merge', merged_bam_file].extend(bamfiles)
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('samtools pileup - step error: %s' % stderr)

    #multiple alignment
    cmd = ['samtools', 'pileup', '-f', reffile, merged_bam_file]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        msg = 'Error:\nstep:%s\n error: %s' % (" ".join(cmd), stderr)
        raise RuntimeError(msg)
    else:
        pileup.write(stdout)
    pileup.close()


if __name__ == '__main__':
    main()