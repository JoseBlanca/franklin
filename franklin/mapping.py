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
from franklin.utils.seqio_utils import seqio
from franklin.sam import (sam2bam, sort_bam_sam, bam2sam, merge_sam,
                          remove_unmapped_reads)
from gff import features_in_gff
import os, shutil

from tempfile import NamedTemporaryFile
from franklin.seq.writers import SamWriter
from franklin.seq.readers import guess_seq_file_format

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

def map_reads_with_tophat(reference_fpath, reads_fpath, out_bam_fpath,
                           parameters):
    'It maps the reads using tophat'
    colorspace = parameters['colorspace']
    threads    = parameters['threads']
    threads    = get_num_threads(threads)

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

def map_reads_with_boina(reference_fpath, reads_fpath, out_bam_fpath,
                         parameters):
    'It maps the reads using first bwa and then tophat'
    # first we are going to map with bwa
    bwa_bam_fhand = NamedTemporaryFile(suffix='.bam')
    bwa_bam_fpath = bwa_bam_fhand.name
    map_reads_with_bwa(reference_fpath, reads_fpath, bwa_bam_fpath, parameters)

    # get unmapped_reads and create new read file
    mapped_bam_fhand        = NamedTemporaryFile(suffix='.bam')
    out_removed_reads_fhand = NamedTemporaryFile(suffix='.sfastq')
    remove_unmapped_reads(open(bwa_bam_fhand.name), mapped_bam_fhand,
                          out_removed_reads_fhand)
    # if the read fileis empty all reads are mapped with bwa
    # it doesn't need to map with tophat
    firstline = open(out_removed_reads_fhand.name).readline()
    if firstline:
        # map with tophat
        tophat_bam_fhand = NamedTemporaryFile(suffix='.bam')
        tophat_bam_fpath = tophat_bam_fhand.name
#        print out_removed_reads_fhand.name
#        raw_input()
        map_reads_with_tophat(reference_fpath, out_removed_reads_fhand.name,
                              tophat_bam_fpath, parameters)

        # merge both bams
        temp_bwa_sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(mapped_bam_fhand.name, temp_bwa_sam_fhand.name, header=True)
        temp_tophat_sam_fhand = NamedTemporaryFile(suffix='.sam')
        bam2sam(tophat_bam_fpath, temp_tophat_sam_fhand.name, header=True)
        temp_out_bam_fhand = NamedTemporaryFile(suffix='.bam')
        merge_sam([open(temp_bwa_sam_fhand.name),
                   open(temp_tophat_sam_fhand.name)],
                   temp_out_bam_fhand, open(reference_fpath))
#        print open(out_bam_fhand.name).read()
#        raw_input()
        tophat_bam_fhand.close()
        temp_tophat_sam_fhand.close()
        temp_bwa_sam_fhand.close()
    else:
        temp_out_bam_fhand = NamedTemporaryFile(suffix='.bam')
        temp_out_bam_fhand.write(open(bwa_bam_fpath).read())
        temp_out_bam_fhand.flush()

    #sort bam
    java_conf = parameters['java_conf']
    sort_bam_sam(temp_out_bam_fhand.name, out_bam_fpath,
                 sort_method='coordinate', java_conf=java_conf)
    #close all temp files
    bwa_bam_fhand.close()
    temp_out_bam_fhand.close()
    out_removed_reads_fhand.close()
    mapped_bam_fhand.close()

def create_gmap_reference(reference_fpath):
    'It creates the reference fpath'
    dir_, name = os.path.split(reference_fpath)
    cmd = ['gmap_setup', '-d', name, '-D', dir_, reference_fpath]
    call(cmd, raise_on_error=True)

    cmd = ['make', '-f', 'Makefile.%s' % name, 'coords']
    call(cmd, raise_on_error=True)

    cmd = ['make', '-f', 'Makefile.%s' % name, 'gmapdb']
    call(cmd, raise_on_error=True)

    cmd = ['make', '-f', 'Makefile.%s' % name, 'install']
    call(cmd, raise_on_error=True)

def map_reads_with_gmap(reference_fpath, reads_fpath, out_bam_fpath,
                        parameters):
    'It maps the reads with gmap'
    threads = parameters['threads']
    reference_dir, reference_name = os.path.split(reference_fpath)
    if not os.path.exists(reference_fpath + '.chromosome'):
        create_gmap_reference(reference_fpath)

    # gmap needs the reads to be fasta and not fastq
    reads_fhand = open(reads_fpath)
    reads_format = guess_seq_file_format(reads_fhand)
    if reads_format != 'fasta':
        out_fhand = NamedTemporaryFile(suffix='.fasta')
        seqio(reads_fhand, out_fhand, 'fasta')
        reads_to_cmap = out_fhand.name
    else:
        reads_to_cmap = reads_fpath

    # create command
    cmd  = ['gmap', '-d', reference_name, '-D', reference_dir, '-f', '4']
    if threads:
        cmd.extend(['-t', str(threads)])
    cmd.append(reads_to_cmap)


    gff3_fhand = NamedTemporaryFile(suffix='.gff3')
    call(cmd, stdout=gff3_fhand, raise_on_error=True)

    out_sam_fhand = NamedTemporaryFile(suffix='.sam')
    gmap_gff_to_sam(open(gff3_fhand.name), open(reference_fpath),
                     open(reads_fpath), out_sam_fhand)
    #print open(out_sam_fhand.name).read()
    sam2bam(out_sam_fhand.name, out_bam_fpath)

    out_sam_fhand.close()
    gff3_fhand.close()

def gmap_gff_to_sam(in_gmap_gff3, reference_fhand, reads_fhand, output_fhand,
                    keep_unmapped=False):
    'It converts the gmap gff3 to sam format'
    samwriter = SamWriter(reference_fhand, reads_fhand, output_fhand,
                          keep_unmapped=keep_unmapped)
    for feature, mapped in features_in_gmap_gff(in_gmap_gff3):
        alignment = gff_feature_to_alignment(feature, mapped)
        samwriter.write(alignment)

def features_in_gmap_gff(in_gmap_gff3):
    '''It returns features in a gmap gff3 and analices if a feature is
    already mapped. The gff must be ordered'''
    old_name     = None
    new_features = []
    for feature in features_in_gff(in_gmap_gff3, version=3):
        feat_name = feature['attributes']['Name']
        if old_name is not None and old_name != feat_name:
            if len(new_features) > 1:
                mapped = False
            else:
                mapped = True
            for new_feature in new_features:
                yield new_feature, mapped
            new_features = []
        new_features.append(feature)
        old_name = feat_name
    else:
        if len(new_features) > 1:
            mapped = False
        else:
            mapped = True
        for new_feature in new_features:
            yield new_feature, mapped


def gff_feature_to_alignment(feature, mapped):
    'converts a gff feature into an alignment struct'
    items       = feature['attributes']['Target'].split()
    query_start = int(items[1])
    strand      = items[3]
    cigar       = _get_cigar_from_feature(feature, query_start)
    alignment = {'reference_name':feature['seqid'],
                 'query_name': feature['attributes']['Name'],
                 'mapped':mapped,
                 'position':feature['start'],
                 'cigar':cigar,
                 'strand':strand}
    return alignment

def _get_cigar_from_feature(feature, query_start):
    'It converts an unclipeed cigar into a clipped cigar alignment'
    cigar = []
    for cigar_item in feature['attributes']['Gap'].split():
        type_ = cigar_item[0]
        num   = cigar_item[1:]
        cigar.append(num + type_)

    if query_start != 1:
        cigar.insert(0, '%dS' % (query_start - 1))

    return ''.join(cigar)

MAPPER_FUNCS = {'bwa': map_reads_with_bwa,
                'tophat':map_reads_with_tophat,
                'gmap':map_reads_with_gmap,
                'boina':map_reads_with_boina}

def map_reads(mapper, reference_fpath, reads_fpath, out_bam_fpath,
              parameters=None):
    'It maps the reads to the reference and returns a bam file'
    if parameters is None:
        parameters = {}
    mapper_func = MAPPER_FUNCS[mapper]
    return mapper_func(reference_fpath, reads_fpath, out_bam_fpath, parameters)
