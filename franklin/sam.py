'''
This library is used to manipulate the sam and bam files. You can merge bams and
add some xtra information to sams

Created on 05/01/2010

@author: peio
'''
from tempfile import NamedTemporaryFile

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

import os, re, tempfile, copy
from itertools import imap
import pysam
from franklin.utils.cmd_utils import call, java_cmd, guess_java_install_dir
from franklin.utils.seqio_utils import seqs_in_file
from franklin.utils.misc_utils import get_num_threads
from franklin.utils.itertools_ import ungroup, classify, take_sample
from franklin.statistics import create_distribution

def get_read_group_info(bam):
    'It returns a dict witht the read group info: platform, lb, etc'
    rg_info = {}
    for read_group in bam.header['RG']:
        name = read_group['ID']
        del read_group['ID']
        rg_info[name] = read_group
    return rg_info

def _guess_picard_path(java_conf=None):
    'It returns the picard path using locate'
    if java_conf and 'picard_path' in java_conf and java_conf['picard_path']:
        return java_conf['picard_path']
    else:
        return guess_java_install_dir('SortSam.jar')

def _guess_gatk_path(java_conf=None):
    'It returns the GATK path using locate'
    if java_conf and 'gatk_path' in java_conf and java_conf['gatk_path']:
        return java_conf['gatk_path']
    else:
        return guess_java_install_dir('GenomeAnalysisTK.jar')

def bam2sam(bam_path, sam_path, header=False):
    '''It converts between bam and sam.'''
    cmd = ['samtools', 'view']
    if header:
        cmd.append('-h')
    cmd.append(bam_path)
    cmd.extend(['-o', sam_path])
    call(cmd, raise_on_error=True)

def sam2bam(sam_path, bam_path):
    'It converts between bam and sam.'
    cmd = ['samtools', 'view', '-bSh', '-o', bam_path, sam_path]
    call(cmd, raise_on_error=True)
    
    #if we try to use this function with a sam with no header it would fail
     #samtools view -bhS -o seqs.bam  seqs.2.sam
    #[samopen] no @SQ lines in the header.
    #[sam_read1] missing header? Abort!
    #in that case we would have to use the -t option

def bamsam_converter(input_fhand, output_fhand, java_conf=None):
    'Converts between sam and bam'
    picard_path = _guess_picard_path(java_conf)
    picard_jar = os.path.join(picard_path, 'SamFormatConverter.jar')
    cmd = java_cmd(java_conf)
    cmd.extend(['-jar', picard_jar, 'INPUT=' + input_fhand,
                'OUTPUT=' + output_fhand])
    call(cmd, raise_on_error=True)

def add_header_and_tags_to_sam(sam_fhand, new_sam_fhand):
    '''It adds tags to each of the reads of the alignments on the bam.
    It creates a new bam in other to avoid errors'''

    # first we write the headers
    sam_fhand.seek(0)
    for line in sam_fhand:
        line = line.strip()
        if not line:
            continue
        if line.startswith('@'):
            new_sam_fhand.write(line + '\n')

    # we add a new header. We take this information from the name of the file
    prefix = os.path.basename(sam_fhand.name)
    prefix = prefix.split('.')[:-1]

    rgid_ = []
    append_to_header = []
    for item in prefix:
        try:
            key, value = item.split('_', 1)
            if key.upper() not in ['LB', 'PL', 'SM']:
                continue
        except ValueError:
            continue
        append_to_header.append('%s:%s' % (key.upper(), value))
        rgid_.append(value)
    rgid = "_".join(rgid_)

    # a proper sam with RG needs a SM key
    if not _check_tag_in_RG('SM', append_to_header):
        append_to_header.append('SM:%s' % rgid)


    new_sam_fhand.write('@RG\tID:%s\t' % rgid)
    new_sam_fhand.write("\t".join(append_to_header) + " \n")

    # now the alignments
    sam_fhand.seek(0)
    for line in sam_fhand:
        line = line.strip()
        if not line:
            continue
        if not line.startswith('@'):
            line = line + "\t" + 'RG:Z:%s' % rgid
            new_sam_fhand.write(line + '\n')
    new_sam_fhand.flush()

def _check_tag_in_RG(tag, tags_to_append):
    'It check if the given tag is in the tags to append'
    for tag_to_append in tags_to_append:
        key = tag_to_append.split(':')[0]
        if key.upper() == tag.upper():
            return True
    return False

def merge_sam(infiles, outfile, reference):
    'It merges a list of sam files'

    #first the reference part of the header
    ref_header = []
    for seq in seqs_in_file(reference):
        name = seq.name
        length = len(seq)
        ref_header.append(['@SQ', 'SN:%s' % name, 'LN:%d' % length])

    #now the read groups
    headers = set()

    for input_ in infiles:
        input_.seek(0)
        for line in input_:
            line = line.strip()
            if line.startswith('@SQ'):
                continue
            elif line.startswith('@'):
                headers.add(tuple(line.split()))
            else:
                break

    #join and write both header parts
    ref_header.extend(headers)
    for header in ref_header:
        outfile.write('\t'.join(header))
        outfile.write('\n')

    #the non header parts
    for input_ in infiles:
        input_.seek(0)
        for line in input_:
            if line.startswith('@'):
                continue
            outfile.write(line)
    outfile.flush()

def sort_bam_sam(in_fhand, out_fhand, sort_method='coordinate',
                 java_conf=None):
    'It sorts a bam file using picard'
    picard_path = _guess_picard_path(java_conf)
    picard_sort_jar = os.path.join(picard_path, 'SortSam.jar')
    java_cmd_ = java_cmd(java_conf)
    java_cmd_.extend(['-jar', picard_sort_jar, 'INPUT=' + in_fhand,
           'OUTPUT=' + out_fhand, 'SORT_ORDER=' + sort_method])
    call(java_cmd_, raise_on_error=True)

def guess_mapped(flag):
    'Giving the flag, guess if the read is mapped or not'
    non_mapped_flags = set()
    if flag == 0:
        mapped = True
    elif flag in non_mapped_flags:  #just a sohortcut
        mapped = False
    else:
        if bin(flag)[-3] == '1':
            mapped = False
            non_mapped_flags.add(flag)
        else:
            mapped = True
    return mapped


def _fix_non_mapped_reads(items):
    'Given a list with sam line it removes MAPQ from non-mapped reads'
    #is the read mapped?
    flag = int(items[1])
    mapped = guess_mapped(flag)
    if not mapped:
        #RNAME = *
        items[2] = '*'
        #POS = 0
        items[3] = '0'
        #MAPQ = 0
        items[4] = '0'
        #CIGAR = *
        items[5] = '*'

def _add_default_quality(items, sanger_read_groups, default_qual):
    'It adds the quality to the sanger reads that do not have it'
    rg_id = filter(lambda x: x.startswith('RG:Z:'), items)
    if rg_id:
        rg_id = rg_id[0].strip().replace('RG:Z:', '')
    else:
        msg = 'Read group is required to add default qualities'
        raise RuntimeError(msg)
    if rg_id not in sanger_read_groups:
        msg = 'Platform for read group %s is not sanger' % rg_id
        raise RuntimeError(msg)
    #here we set the default qual
    items[10] = default_qual * len(items[9])

def standardize_sam(in_fhand, out_fhand, default_sanger_quality=None,
                   add_def_qual=False, only_std_char=True,
                   fix_non_mapped=True):
    '''It adds the default qualities to the reads that do not have one and
    it makes sure that only ACTGN are used in the sam file reads.

    The quality should be given as a phred integer
    '''
    in_fhand.seek(0)
    sep = '\t'

    sanger_read_groups = set()
    if default_sanger_quality:
        default_qual = chr(int(default_sanger_quality) + 33)
    else:
        default_qual = None
    regex = re.compile(r'[^ACTGN]')
    for line in in_fhand:
        #it we're adding qualities we have to keep the read group info
        if add_def_qual and line[:3] == '@RG':
            items = line.split(sep)
            rg_id = None
            platform = None
            for item in items:
                if item.startswith('PL:'):
                    platform = item[3:].strip().lower()
                elif item.startswith('ID:'):
                    rg_id = item[3:].strip()
            if platform == 'sanger':
                sanger_read_groups.add(rg_id)
        elif line[0] == '@':
            pass
        else:
            items = line.split(sep)
            #quality replaced
            do_qual = add_def_qual and items[10] == '*'
            if do_qual or only_std_char:
                #we add default qualities to the reads with no quality
                if do_qual: #empty quality
                    _add_default_quality(items, sanger_read_groups,
                                         default_qual)
                #we change all seq characters to ACTGN
                if only_std_char:
                    text = items[9].upper()
                    items[9] = regex.sub('N', text)
                #we remove MAPQ from the non-mapped reads
                if fix_non_mapped:
                    _fix_non_mapped_reads(items)
                line = sep.join(items)
        out_fhand.write(line)

def create_bam_index(bam_fpath):
    'It creates an index of the bam if it does not exist'
    index_fpath = bam_fpath + '.bai'

    if os.path.exists(index_fpath):
        return
    cmd = ['samtools', 'index', bam_fpath]
    call(cmd, raise_on_error=True)

def create_picard_dict(reference_fpath, java_conf=None):
    'It creates a picard dict if if it does not exist'
    dict_path = os.path.splitext(reference_fpath)[0] + '.dict'
    if os.path.exists(dict_path):
        return
    picard_path = _guess_picard_path(java_conf)
    picard_jar = os.path.join(picard_path, 'CreateSequenceDictionary.jar')
    cmd = ['java', '-jar', picard_jar,
           'R=%s' % reference_fpath,
           'O=%s' % dict_path]
    call(cmd, raise_on_error=True)

def create_sam_reference_index(reference_fpath):
    'It creates a sam index for a reference sequence file'
    index_fpath = reference_fpath + '.fai'
    if os.path.exists(index_fpath):
        return
    cmd = ['samtools', 'faidx', reference_fpath]
    call(cmd, raise_on_error=True)


def realign_bam(bam_fpath, reference_fpath, out_bam_fpath, java_conf=None,
                threads=False):
    'It realigns the bam using GATK Local realignment around indels'
    #reference sam index
    create_sam_reference_index(reference_fpath)

    #reference picard dict
    create_picard_dict(reference_fpath, java_conf=java_conf)

    #bam index
    create_bam_index(bam_fpath)

    #the intervals to realign
    gatk_path = _guess_gatk_path(java_conf)
    gatk_jar = os.path.join(gatk_path, 'GenomeAnalysisTK.jar')
    intervals_fhand = tempfile.NamedTemporaryFile(prefix='intervals',
                                                  suffix='.txt')
    cmd = java_cmd(java_conf=java_conf)
    cmd.extend(['-jar', gatk_jar, '-T', 'RealignerTargetCreator',
           '-I', bam_fpath, '-R', reference_fpath, '-o', intervals_fhand.name])
    #according to GATK this is experimental, so it might be a good idea to
    #do it in just one thread
    parallel = True
    if parallel:
        cmd.extend(['--numthreads', str(get_num_threads(threads))])
    call(cmd, raise_on_error=True)

    #the realignment itself
    cmd = java_cmd(java_conf=java_conf)
    cmd.extend(['-Djava.io.tmpdir=%s' % tempfile.gettempdir(),
           '-jar', gatk_jar, '-I', bam_fpath, '-R', reference_fpath,
           '-T', 'IndelRealigner', '-targetIntervals', intervals_fhand.name,
           '--output', out_bam_fpath])
    if parallel:
        cmd.extend(['--numthreads', str(get_num_threads(threads))])
    call(cmd, raise_on_error=True)

def _calculate_bam_coverage_data(bam):
    '''This function gets data to make stats of a sam file.

    It extracts per column coverage data from a bam
    '''
    # TODO, needs sampling
    for column in bam.pileup():
        reads_per_colum = {}
        for pileup_read in column.pileups:
            aligned_read = pileup_read.alignment
            read_group = aligned_read.opt('RG')
            if read_group not in reads_per_colum:
                reads_per_colum[read_group] = 0
            reads_per_colum[read_group] += 1

        yield reads_per_colum

def _get_bam_coverage(bam, rgs, grouping):
    'It serializes the coverage_data'
    coverages = _calculate_bam_coverage_data(bam)
    coverages = ungroup(coverages, lambda x: x.items())
    coverages = imap(lambda x: (rgs[x[0]][grouping], x[1]), coverages)
    coverages = classify(coverages, lambda x: x[0])
    new_coverages = {}
    for group, values in coverages.items():
        values = imap(lambda x: x[1], values)
        new_coverages[group] = values

    return new_coverages

def _generate_bam_mapping_quality_data(bam):
    '''This function get data to make stats of a sam file.

    It extracts mapping quality per read
    '''
    for aligned_read in pysam.IteratorRowAll(bam):
        read_mapping_qual = aligned_read.mapq
        read_group = aligned_read.opt('RG')
        yield(read_group, read_mapping_qual)

def _get_bam_mapping_quality(bam, rgs, grouping):
    'It serializes the mapq data'
    mapqs = _generate_bam_mapping_quality_data(bam)
    mapqs = imap(lambda x: (rgs[x[0]][grouping], x[1]), mapqs)
    mapqs = classify(mapqs, lambda x: x[0])
    new_mapqs = {}
    for group, values in mapqs.items():
        values = imap(lambda x: x[1], values)
        new_mapqs[group] = values

    return new_mapqs

def bam_distribs(bam_fhand, kind, basename=None, range_=None, grouping=None,
                 sample_size=None):
    '''It makes the bam coverage distribution.

    It can make the distribution taking into account any of the readgroup items:
    platform, sample and library
    '''
    value_calculator = {'coverage':_get_bam_coverage,
                       'mapq':_get_bam_mapping_quality}

    coverage_labels = {'title': "Coverage for %s %s",
                       'xlabel': 'Coverage' ,
                       'ylabel': 'Num. of positions'
                       }
    mapping_labels = {'title': "Mapping qualities for %s %s",
                       'xlabel': "mapping quality" ,
                       'ylabel': 'Num. of reads'
                       }
    plot_labels = {'coverage': coverage_labels,
                   'mapq':mapping_labels}

    if sample_size is not None:
        sampled_bam_fhand = NamedTemporaryFile(suffix='.bam')
        sample_bam(bam_fhand, sampled_bam_fhand, sample_size)
        sample_fpath = sampled_bam_fhand.name
    else:
        sample_fpath = bam_fhand.name

    create_bam_index(bam_fpath=sample_fpath)
    bam = pysam.Samfile(sample_fpath, 'rb')
    rgs = get_read_group_info(bam)
    if grouping is None:

        platforms = set([rg['PL'] for rg in rgs.values()])
        if len(platforms) > 1:
            grouping = 'PL'
        else:
            grouping = 'SM'
    item_values = value_calculator[kind](bam, rgs, grouping)
    results = {}
    for group_name, values in item_values.items():
        if basename is None:
            distrib_fhand = None
            plot_fhand = None
        else:
            distrib_fhand = open('%s.%s_%s.dat' % (basename, kind, group_name),
                                 'w')
            plot_fhand = open('%s.%s_%s.png' % (basename, kind, group_name),
                              'w')
        if grouping == 'PL':
            grouping = 'platform'
        elif grouping == 'SM':
            grouping = 'sample'
        labels = copy.deepcopy(plot_labels[kind])
        labels['title'] = labels['title'] % (grouping, group_name)

        range_ = range_

        distrib = create_distribution(values, labels, distrib_fhand=distrib_fhand,
                                   plot_fhand=plot_fhand, range_=range_)

        results[(grouping, group_name)] = distrib['distrib']
    return results

def sample_bam(bam_fhand, out_bam_fhand, sample_size):
    'It takes a sample from a bam'
    sam_fhand = NamedTemporaryFile(suffix='.sam')
    final_sam = NamedTemporaryFile(suffix='.sam')
    bam2sam(bam_fhand.name, sam_fhand.name, header=True)

    # First get header
    for line in open(sam_fhand.name):
        if line[0] == '@':
            final_sam.write(line)
        else:
            break
    sam_body = take_sample(_reads_in_sam(sam_fhand), sample_size=sample_size)

    for line in sam_body:
        final_sam.write(line)
    final_sam.flush()
    sam2bam(final_sam.name, out_bam_fhand.name)

def _reads_in_sam(sam_fhand):
    'It yields a read mapping in each iteration'
    for line in open(sam_fhand.name):
        if line[0] == '@':
            continue
        yield line
