'''
Created on 05/02/2010

@author: peio
'''
from franklin.utils.cmd_utils import call
from franklin.utils.misc_utils import NamedTemporaryDir
import os

def create_bwa_reference(reference_fpath):
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
    call(cmd, raise_on_error=True)

def map_reads_with_bwa(reference_fpath, reads_fpath, bam_fpath,
                       parameters):
    'It maps the reads to the reference using bwa and returns a bam file'
    #the reference should have an index
    bwt_fpath = os.path.basename(reference_fpath) + '.bwt'
    if not os.path.exists(bwt_fpath):
        create_bwa_reference(reference_fpath)
    reads_length = parameters['reads_length']

    temp_dir = NamedTemporaryDir()
    if reads_length == 'short':
        cmd = ['bwa', 'aln', reference_fpath, reads_fpath]
        sai_fhand = open(os.path.join(temp_dir.name, 'output.sai'), 'wb')
        call(cmd, stdout=sai_fhand, raise_on_error=True)

        cmd = ['bwa', 'samse', reference_fpath, sai_fhand.name, reads_fpath]
        ali_fhand = open(os.path.join(temp_dir.name, 'output.ali'), 'w')
        call(cmd, stdout=ali_fhand, raise_on_error=True)
    elif reads_length == 'long':
        cmd = ['bwa', 'dbwtsw', reference_fpath, reads_fpath]
        ali_fhand = open(os.path.join(temp_dir.name, 'output.ali'), 'w')
        call(cmd, stdout=ali_fhand, raise_on_error=True)
    else:
        raise ValueError('Reads length: short or long')
    # From sam to Bam
    unsorted_bam = os.path.join(temp_dir.name, 'bam_file.bam')
    cmd = ['samtools', 'view' , '-bt', reference_fpath, '-o', unsorted_bam,
           ali_fhand.name]
    call(cmd, raise_on_error=True)
    # sort bam file
    bam_basename = os.path.splitext(bam_fpath)[0]
    cmd = ['samtools', 'sort', unsorted_bam, bam_basename]
    call(cmd, raise_on_error=True)

    temp_dir.close()


MAPPER_FUNCS = {'bwa': map_reads_with_bwa}

def map_reads(mapper, reference_fpath, reads_fpath, out_bam_fpath,
              parameters=None):
    'It maps the reads to the reference and returns a bam file'
    if parameters is None:
        parameters = {}
    mapper_func = MAPPER_FUNCS[mapper]
    return mapper_func(reference_fpath, reads_fpath, out_bam_fpath, parameters)
