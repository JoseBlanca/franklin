'''
Created on 05/02/2010

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
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

from franklin.utils.cmd_utils import call
from franklin.utils.misc_utils import NamedTemporaryDir, get_num_threads
from franklin.sam import sam2bam, sort_bam_sam
import os, shutil

def create_bwa_reference(reference_fpath, color=False):
    'It creates the bwa index for the given reference'
    #how many sequences do we have?
    n_seqs = 0
    for line in open(reference_fpath):
        if line[0] == '>':
            n_seqs += 1
    if n_seqs > 10000:
        algorithm = 'bwtsw'
    else:
        algorithm = 'is'

    cmd = ['bwa', 'index', '-a', algorithm, reference_fpath]
    if color:
        cmd.append('-c')

    call(cmd, raise_on_error=True)


def map_reads_with_bwa(reference_fpath, reads_fpath, bam_fpath,
                       parameters, threads=False, java_conf=None):
    'It maps the reads to the reference using bwa and returns a bam file'
    colorspace   = parameters['colorspace']
    reads_length = parameters['reads_length']
    threads = get_num_threads(threads)
    #the reference should have an index
    bwt_fpath = reference_fpath + '.bwt'
    if not os.path.exists(bwt_fpath):
        create_bwa_reference(reference_fpath, color=colorspace)

    temp_dir = NamedTemporaryDir()
    output_ali = 'output.ali'
    bam_file_bam = 'bam_file.bam'
    output_sai = 'output.sai'
    if reads_length == 'short':
        cmd = ['bwa', 'aln', reference_fpath, reads_fpath,
               '-t', str(threads)]
        if colorspace:
            cmd.append('-c')
        sai_fhand = open(os.path.join(temp_dir.name, output_sai), 'wb')
        call(cmd, stdout=sai_fhand, raise_on_error=True)

        cmd = ['bwa', 'samse', reference_fpath, sai_fhand.name, reads_fpath]
        ali_fhand = open(os.path.join(temp_dir.name, output_ali), 'w')
        call(cmd, stdout=ali_fhand, raise_on_error=True)
    elif reads_length == 'long':
        cmd = ['bwa', 'dbwtsw', reference_fpath, reads_fpath,
               '-t', str(threads)]
        ali_fhand = open(os.path.join(temp_dir.name, output_ali), 'w')
        call(cmd, stdout=ali_fhand, raise_on_error=True)
    else:
        raise ValueError('Reads length: short or long')
    # From sam to Bam
    unsorted_bam = os.path.join(temp_dir.name, bam_file_bam)
    sam2bam(ali_fhand.name, unsorted_bam)
    # sort bam file
    sort_bam_sam(unsorted_bam, bam_fpath, sort_method='coordinate',
                 java_conf=java_conf, strict_validation=False)
    temp_dir.close()

def create_bowtie_reference(reference_fpath):
    'It creates the bowtie index used by bowtie and tophat'
    bowtie_index = os.path.splitext(reference_fpath)[0]
    cmd = ['bowtie-build', reference_fpath, bowtie_index]
    call(cmd, raise_on_error=True)

def map_reads_with_tophat(reference_fpath, reads_fpath, out_bam_fpath,
                           parameters, threads=None, java_conf=None):
    'It maps the reads using tophat'
    threads      = get_num_threads(threads)
    colorspace   = parameters['colorspace']

    # create needed index
    bowtie_index = os.path.splitext(reference_fpath)[0]
    if not os.path.exists(bowtie_index + '.1.ebwt'):
        create_bowtie_reference(reference_fpath)

    # run tophat
    temp_dir = NamedTemporaryDir()
    cmd      = ['tophat', '-o', temp_dir.name]
    if colorspace:
        cmd.append('-C')
    if threads:
        cmd.extend(['-p', str(threads)])

    cmd.extend([bowtie_index, reads_fpath])
    call(cmd, raise_on_error=True)

    # copy the output
    shutil.copy(os.path.join(temp_dir.name, 'accepted_hits.bam'), out_bam_fpath)
    temp_dir.close()

MAPPER_FUNCS = {'bwa': map_reads_with_bwa,
                'tophat':map_reads_with_tophat}

def map_reads(mapper, reference_fpath, reads_fpath, out_bam_fpath,
              parameters=None, threads=False, java_conf=None):
    'It maps the reads to the reference and returns a bam file'
    if parameters is None:
        parameters = {}
    mapper_func = MAPPER_FUNCS[mapper]
    return mapper_func(reference_fpath, reads_fpath, out_bam_fpath, parameters,
                       threads=threads, java_conf=java_conf)
