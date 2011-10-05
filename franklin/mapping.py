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

import os, shutil
from tempfile import NamedTemporaryFile

from franklin.utils.cmd_utils import call, get_external_bin_dir
from franklin.utils.misc_utils import NamedTemporaryDir, get_num_threads
from franklin.sam import (sam2bam, sort_bam_sam)

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


def map_reads_with_bwa(reference_fpath, reads_fpath, bam_fpath, parameters):
    'It maps the reads to the reference using bwa and returns a bam file'
    colorspace   = parameters['colorspace']
    reads_length = parameters['reads_length']
    threads      = parameters['threads']
    java_conf    = parameters['java_conf']
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

def create_bowtie_reference(reference_fpath, color=False):
    'It creates the bowtie index used by bowtie and tophat'
    bowtie_index = os.path.splitext(reference_fpath)[0]
    cmd = ['bowtie-build']
    if color:
        cmd.append('-C')

    cmd.extend([reference_fpath, bowtie_index])
    call(cmd, raise_on_error=True)

def _modify_gmap_makefile(makefile):
    'It modifies the gmap makefile to use proper binaries'
    fhand = NamedTemporaryFile(delete=False)
    bin_dir = get_external_bin_dir()
    for line in open(makefile):
        if line.startswith('FA_COORDS'):
            line = 'FA_COORDS = %s/fa_coords\n' % bin_dir
        if line.startswith('MD_COORDS'):
            line = 'MD_COORDS = %s/md_coords\n' % bin_dir
        if line.startswith('GMAP_PROCESS'):
            line = 'GMAP_PROCESS = %s/gmap_process\n' % bin_dir
        if line.startswith('GMAPINDEX'):
            line = 'GMAPINDEX = %s/gmapindex\n' % bin_dir
        if line.startswith('PMAPINDEX'):
            line = 'PMAPINDEX = %s/pmapindex\n' % bin_dir
        fhand.write(line)
    fhand.flush()
    shutil.move(fhand.name, makefile)

def create_gmap_reference_old(reference_fpath):
    'It creates the reference fpath'
    dir_, name = os.path.split(reference_fpath)
    if not dir_:
        dir_ = '.'
    makefile_fpath = os.path.join(dir_, 'Makefile.%s' % name)
    #the gmap_setup command would not accept a file, to avoid thread conflicts
    #we first create a name with a random name and them we move it where it
    #belongs
    temp_makefile = NamedTemporaryFile(delete=False)
    cmd = ['gmap_setup', '-d', name, '-D', dir_, '-o', temp_makefile.name,
           reference_fpath]
    try:
        call(cmd, raise_on_error=True)
    except OSError:
        raise OSError('Gmap mapper is not installed or not in the path')
    shutil.move(temp_makefile.name, makefile_fpath)

    _modify_gmap_makefile(makefile_fpath)

    cmd = ['make', '-f', makefile_fpath , 'coords']
    call(cmd, raise_on_error=True, add_ext_dir=False)

    cmd = ['make', '-f', makefile_fpath, 'gmapdb']
    call(cmd, raise_on_error=True, add_ext_dir=False)

    cmd = ['make', '-f', makefile_fpath, 'install']
    call(cmd, raise_on_error=True, add_ext_dir=False)

    #we remove the makefile and an extra file with some instructions
    os.remove(makefile_fpath)
    coords_fpath = 'coords.%s' % name
    install_fpath = 'INSTALL.%s' % name
    for fpath in (coords_fpath, install_fpath):
        if os.path.exists(fpath):
            os.remove(fpath)

def create_gmap_reference(reference_dir, reference_path, reference_name,
                          parameters=None):
    'It creates the reference fpath'

    cmd = ['gmap_build',  '-B',  get_external_bin_dir(), '-D',  reference_dir,
           '-d', reference_name]
    if parameters and 'kmer' in parameters:
        cmd.extend(['-k', str(parameters['kmer'])])
    cmd.append(reference_path)
    call(cmd,raise_on_error=True)

    fpath = '%s.coords' % reference_name
    if os.path.exists(fpath):
        os.remove(fpath)

def map_reads_with_gmap(reference_fpath, reads_fpath, out_bam_fpath,
                        parameters):
    'It maps the reads with gmap'
    threads = parameters['threads']
    reference_dir, reference_file_name = os.path.split(reference_fpath)
    reference_name = reference_file_name.split('.')[0]
    if not reference_dir:
        reference_dir = '.'

    if not os.path.exists(os.path.join(reference_dir, reference_name,
                                       reference_name + '.chromosome')):
        create_gmap_reference(reference_dir, reference_fpath, reference_name,
                              parameters)

    cmd  = ['gmap', '-d', reference_name, '-D', reference_dir, '-f', 'samse']
    # this gmap options doesn' detect deletions close to introns
    cmd.append('--canonical-mode=0')
    if threads:
        cmd.extend(['-t', str(threads)])
    cmd.append(reads_fpath)
    out_sam_fhand = NamedTemporaryFile(suffix='.sam')
    call(cmd, stdout=out_sam_fhand, raise_on_error=True)
    sam2bam(out_sam_fhand.name, out_bam_fpath)
    out_sam_fhand.close()


MAPPER_FUNCS = {'bwa': map_reads_with_bwa,
                'gmap':map_reads_with_gmap}

def map_reads(mapper, reference_fpath, reads_fpath, out_bam_fpath,
              parameters=None):
    'It maps the reads to the reference and returns a bam file'
    if parameters is None:
        parameters = {}
    mapper_func = MAPPER_FUNCS[mapper]
    return mapper_func(reference_fpath, reads_fpath, out_bam_fpath, parameters)
