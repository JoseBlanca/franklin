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
It runs bwa.

It needs tha the database is indexed. It can run any type of sequence

Created on 08/10/2009

@author: peio
'''
from optparse import OptionParser
import os
from biolib.biolib_cmd_utils import NamedTemporaryDir, call

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i sam_pileup, -o req_pos -p pipeline')
    parser.add_option('-s', '--seqs', dest='seqfile', help='Sequences file')
    parser.add_option('-f', '--type', dest='seq_type', default='short',
                      help='Sequences type: short or long')
    parser.add_option('-r', '--references', dest='reffile',
                      help='reference seq file')
    parser.add_option('-i', '--refdb', dest='refdb', help='reference index')
    parser.add_option('-b', '--bamfile', dest='bamfile', help='bam output file',
                      default='output')

    return parser
def set_parameters():
    '''It sets the parameters for the script.'''

    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.seqfile is None:
        parser.error('seqfile is mandatory')
    else:
        seqfile = options.seqfile

    if options.reffile is None:
        parser.error('file is mandatory')
    else:
        reffile = options.reffile

    seq_type = options.seq_type
    refdb    = options.refdb
    bamfile  = options.bamfile



    return seqfile, reffile, seq_type, refdb, bamfile

def main():
    "The real script"
    seqfile, reffile, seq_type, refdb, bamfile = set_parameters()

    # we need that the index is created
    if not os.path.exists(refdb + '.bwt'):
        raise OSError('Reference database index not created: %s' % refdb)

    # TEMPORARY DIR NAME
    temp_dir = NamedTemporaryDir()
    dir_name = temp_dir.name

    if seq_type == 'short':
        cmd = ['bwa', 'aln', refdb, seqfile]
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('bwa aln -step error: %s' % stderr)
        sai_fhand = open(os.path.join(dir_name, 'output.sai'), 'wb')
        sai_fhand.write(stdout)
        sai_fhand.close()

        cmd = ['bwa', 'samse', refdb, sai_fhand.name, seqfile]
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('bwa samse - step error: %s' % stderr)
        ali_fhand = open(os.path.join(dir_name, 'output.ali'), 'w')
        ali_fhand.write(stdout)
        ali_fhand.close()

    elif seq_type == 'long':
        #align sanger
        cmd = ['bwa', 'dbwtsw', refdb, seqfile]
        stdout, stderr, retcode = call(cmd)
        if retcode:
            raise RuntimeError('bwa dbwtsw - step error: %s\n%s' % (stderr, " ".join(cmd)))
        ali_fhand = open(os.path.join(dir_name, 'output.ali'), 'w')
        ali_fhand.write(stdout)
        ali_fhand.close()
    else:
        raise ValueError('Seq type: short or long')

    #from sam import to bam
    cmd = ['samtools', 'view' , '-bt', reffile, '-o',
           os.path.join(dir_name, 'bam_file.bam'), ali_fhand.name]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('samtools view - step error: %s' % stderr)

    #sort the bam
    cmd = ['samtools', 'sort', os.path.join(dir_name, 'bam_file.bam'), bamfile]
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError('samtools sort - step error: %s' % stderr)

if __name__ == '__main__':
    main()