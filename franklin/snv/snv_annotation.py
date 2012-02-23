'''
Created on 16/02/2010

@author: jose
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

from __future__ import division
from operator import attrgetter
from collections import defaultdict
from copy import copy
import math, itertools

try:
    import pysam
except ImportError:
    pass

from Bio.SeqFeature import FeatureLocation
from Bio.Restriction import Analysis, CommOnly, RestrictionBatch
from franklin.seq.seqs import SeqFeature, get_seq_name
from franklin.utils.misc_utils import get_fhand
from franklin.sam import create_bam_index, get_read_group_info

DEFAUL_MIN_NUM_READS_PER_ALLELE = 2
DEFAULT_PLOIDY = 2

DELETION_ALLELE = '-'
UNKNOWN_NUCLEOTYDE = 'N'
NO_NUCLEOTYDE = ' '
N_ALLLELES = ('n', '?')

SNP = 0
INSERTION = 1
DELETION = 2
INVARIANT = 3
INDEL = 4
COMPLEX = 5
TRANSITION = 6
TRANSVERSION = 7
UNKNOWN = 8

SNV_TYPES = {SNP:'SNP', INSERTION:'insertion', DELETION:'deletion',
             INVARIANT:'invariant', INDEL:'indel', COMPLEX:'complex',
             TRANSITION:'transition', TRANSVERSION:'transversion',
             UNKNOWN:'unknown'}

COMMON_ENZYMES = ['EcoRI', 'SmaI', 'BamHI', 'AluI', 'BglII',
                  'SalI', 'BglI', 'ClaI', 'TaqI',
                  'PstI', 'PvuII', 'HindIII', 'EcoRV',
                  'HaeIII', 'KpnI', 'ScaI',
                  'HinfI', 'DraI', 'ApaI', 'BstEII', 'ZraI', 'BanI', 'Asp718I']
UNKNOWN_RG = 'unknown'

MATCH = 0
INSERTION = 1
DELETION = 2
SKIP = 3
SOFT_CLIP = 4
HARD_CLIP = 5
PADDING = 6

IN_FIRST_POS = 1
IN_MIDDLE_POS = 2
IN_LAST_POS = 3
IN_FIRST_AND_LAST = 4


CIGAR_DECODE = 'MIDNSHP'


def _get_raw_allele_from_read(aligned_read, index):
    'It returns allele, quality, is_reverse'
    allele = aligned_read.seq[index].upper()
    if aligned_read.qual:
        try:
            qual = _quality_to_phred(aligned_read.qual[index])
        except ZeroDivisionError:
            raise RuntimeError('mi error')
    else:
        qual = None
    return allele, qual

SEGMENTS_CACHE = {}
MAX_CACHED_SEGMENTS = 10000

def _get_cigar_segments_from_aligned_read(aligned_read):
    'It gets the cigar from an aligened_read of pysam'
    begin_pos_read_in_ref = aligned_read.pos
    read_len = len(aligned_read.seq)
    cigar = aligned_read.cigar
    return _get_segments_from_cigar(begin_pos_read_in_ref, cigar, read_len)

def _get_segments_from_cigar(begin_pos_read_in_ref, cigar, read_len):
    '''It returns two lists (reference and read) in which the firsts nucleotides
     of the different cigar categories are given.

     67890  12345
     ATCGA--GATCG
     atcGATCG--CG
        01234  56

    CIGAR = 3H2M2I1M2D2M

    ref_segments = [9, None, 11, 12, 14]
    read_segments = [0, 2, 4, None, 5]

    It also returns the limits of the aligned reference and the limits of the
    read.

    It also returns a list with the cigar category for each segment.
    The position in the reference is only used for the cache
    '''

    cigar = tuple(cigar)
    global SEGMENTS_CACHE
    cache_key = read_len, begin_pos_read_in_ref, cigar
    if cache_key in SEGMENTS_CACHE:
        return SEGMENTS_CACHE[cache_key]['result']
    #We ignore hard clipped nucleotides ('H')
    cigar_elements = []
    for element in range(len(cigar)):
        if cigar[element][0] != HARD_CLIP:
            cigar_elements.append(cigar[element])

    ref_segments = []
    read_segments = []
    ref_pos = begin_pos_read_in_ref
    read_start = 0
    read_end = 0
    segment_type = []
    segment_lens = []

    read_pos = read_start

    for index, element in enumerate(cigar_elements):
        if element[0] == SOFT_CLIP:
            #Is the soft clip at beginning?
            if index == 0:
                read_pos += element[1]
                read_start += element[1]
                read_end = read_start
            elif index == len(cigar_elements) - 1:
                continue
            else:
                msg = 'Soft clips in the middle of the read are not supported'
                raise RuntimeError(msg)
        elif element[0] == MATCH:
            segment_type.append(MATCH)
            ref_segments.append(ref_pos)
            read_segments.append(read_pos)
            segment_lens.append(element[1])
            ref_pos += element[1]
            read_pos += element[1]
            read_end += element[1]
        elif element[0] == INSERTION:
            segment_type.append(INSERTION)
            ref_segments.append(None)
            read_segments.append(read_pos)
            segment_lens.append(element[1])
            read_pos += element[1]
            read_end += element[1]
        elif element[0] == DELETION or element[0] == SKIP:
            read_segments.append(None)
            ref_segments.append(ref_pos)
            segment_lens.append(element[1])
            #We differentiate between DELETION or SKIP
            if element[0] == DELETION:
                segment_type.append(DELETION)
            if element[0] == SKIP:
                segment_type.append(SKIP)

            ref_pos += element[1]

    ref_end = ref_pos - 1
    ref_start = ref_segments[0]
    # Somethimes the first segment is an insertion and the reference por the
    # start position is in the next segment
    if ref_segments[0] is not None:
        ref_start = ref_segments[0]
    else:
        ref_start = ref_segments[1]

    ref_limits = [ref_start, ref_end]

    read_end = read_end - 1
    read_limits = [read_start, read_end]

    result = (ref_segments, read_segments, sorted(ref_limits),
              sorted(read_limits), segment_type, segment_lens)
    if ref_pos is not None:
        #store cache
        SEGMENTS_CACHE[cache_key] = {'result':result, 'ref_pos':ref_pos}
        #clean cache
        if len(SEGMENTS_CACHE) > MAX_CACHED_SEGMENTS:
            SEGMENTS_CACHE = {}
    return result


def _locate_segment(ref_pos, ref_segments, segment_lens, ref_limits):
    'It locates a read position in the segments'

    for segment_index, segment_begin_pos in enumerate(ref_segments):
        if segment_begin_pos is None:
            continue
        end_segment_pos = segment_begin_pos + segment_lens[segment_index] - 1
        if ref_pos < segment_begin_pos:
            #we're before the first segment, no read here
            return None
        elif ref_pos > end_segment_pos:
            continue
        elif ref_pos == segment_begin_pos and segment_begin_pos == end_segment_pos:
            return segment_index, IN_FIRST_AND_LAST
        elif ref_pos == segment_begin_pos:
            return segment_index, IN_FIRST_POS
        elif ref_pos == end_segment_pos:
            return segment_index, IN_LAST_POS
        elif segment_begin_pos < ref_pos < end_segment_pos:
            return segment_index, IN_MIDDLE_POS
        else:
            raise RuntimeError('We should not be here, fix me')
    #we're outside any segment
    return None

def _get_insertion(segment_index, segment_type, read_pos, aligned_read,
                   segment_lens):
    #TODO explain function
    allele = None
    kind = None
    qual = None

    if segment_index == len(segment_type) - 1:
        #we're in the last segment
        next_segment_pos = None
    else:
        next_segment_pos = segment_index + 1


    if (next_segment_pos is not None and
        segment_type[next_segment_pos] == INSERTION):
        indel_length = segment_lens[next_segment_pos]
        if indel_length == 0:
            msg = "An insertion can't be of 0 length\n"
            msg += 'next_segment_pos ' + str(next_segment_pos)
            msg += '\nsegment_index ' + str(segment_index)
            msg += '\nsegment_type ' + str(segment_type)
            msg += '\nread_pos ' + str(read_pos)
            msg += '\naligned_read ' + str(aligned_read)
            raise ValueError(msg)
        start = read_pos
        end = start + indel_length
        allele, qual = _get_raw_allele_from_read(aligned_read,
                                                 slice(start, end))
        kind = INSERTION

    return allele, kind, qual

def _from_ref_to_read_pos(segment_type, ref_segment_pos, read_segment_pos,
                         ref_pos):
    'Given the segment positions it calculates the read_pos'
    if segment_type != MATCH:
        return None
    read_pos = read_segment_pos + (ref_pos - ref_segment_pos)
    return read_pos

def _read_pos_around_del(ref_segments, read_segments, segment_types,
                        segment_index, ref_pos):
    'It returns the read positions around a deletion'
    if segment_types[segment_index - 1] == MATCH:
        read_pos1 = _from_ref_to_read_pos(segment_types[segment_index - 1],
                                         ref_segments[segment_index - 1],
                                         read_segments[segment_index - 1],
                                         ref_segments[segment_index] - 1)
        read_pos2 = read_pos1 + 1
    elif segment_types[segment_index + 1] == MATCH:
        read_pos2 = _from_ref_to_read_pos(segment_types[segment_index + 1],
                                         ref_segments[segment_index + 1],
                                         read_segments[segment_index + 1],
                                         ref_segments[segment_index + 1])
        read_pos1 = read_pos2 - 1
    else:
        msg = 'A deletion is surrounded by two segments that are not matches'
        raise RuntimeError(msg)
    return read_pos1, read_pos2

def _chop_alignment(alignment):
    'it chops the alignments in lists'
    ref, read = alignment
    listed_ref = []
    listed_read = []
    inserts = ''
    read_inserts = ''
    for index in range(len(ref)):
        if ref[index] == '-':
            inserts += '-'
            read_inserts += read[index]
        else:
            if inserts != '':
                listed_ref.append(inserts)
                listed_read.append(read_inserts)
                inserts = ''
                read_inserts = ''
            listed_ref.append(ref[index])
            listed_read.append(read[index])
    else:
        if inserts != '':
            listed_ref.append(inserts)
            listed_read.append(read_inserts)
    return listed_ref, listed_read

def _prepare_alignments(alignments):
    'It prepares the alignements for joining'
    prepared_alignments = []
    for name, alignment in alignments.items():
        ref, read = _chop_alignment(alignment)
        assert len(ref) == len(read)
        alignment = {'reference':ref,
                     'reads':{name:read}}
        prepared_alignments.append(alignment)
    return prepared_alignments

def add_collumn(alignment, nalig, index, diff_length=0):
    'It adds a column to the given alignment'
    for name, read in alignment['reads'].items():
        if name not in nalig['reads']:
            nalig['reads'][name] = []
        try:
            to_add = read[index]
        except IndexError:
            to_add = diff_length * '-'
        nalig['reads'][name].append(to_add)

def _join_alignments(alignment1, alignment2, snv_types_per_read):
    'It joins two alignemnts and makes a new alignment'
    def _insert_in_seq(list_seq):
        for element in list_seq:
            if '-' in element:
                return True
        return False
    if len(alignment1['reference']) > len(alignment2['reference']):
        align1 = alignment1
        align2 = alignment2
    else:
        align1 = alignment2
        align2 = alignment1

    ref1 = align1['reference']
    ref2 = align2['reference']

    if (not _insert_in_seq(ref1) and not _insert_in_seq(ref2) and
        len(ref1) == len(ref2)):
        reads = {}
        for name, read in align1['reads'].items():
            reads[name] = read
        for name, read in align2['reads'].items():
            reads[name] = read
        return {'reference':ref1,
                'reads':reads}

    nalig = {'reference':[], 'reads':{}}
    index1_delta = 0
    index2_delta = 0
    for index in range(len(ref1)):

        index1 = index - index1_delta
        index2 = index - index2_delta
        ref1_item = ref1[index1]
        try:
            ref2_item = ref2[index2]
        except IndexError:
            ref2_item = ref1_item

        inser_1 = ref1_item.count('-')
        inser_2 = ref2_item.count('-')
#        print ref1_item, ref2_item, inser_1, inser_2
        if inser_1 == inser_2:

            nalig['reference'].append(ref1_item)
            add_collumn(align1, nalig, index1, inser_1)
            add_collumn(align2, nalig, index2, inser_2)

        elif inser_1 > inser_2:
            inser_diff = inser_1 - inser_2
            nalig['reference'].append(ref1_item)
            add_collumn(align1, nalig, index1)
            just_once = True
            for name, read in align2['reads'].items():
                if name not in nalig['reads']:
                    nalig['reads'][name] = []
                if inser_2 > 0:
                    nalig['reads'][name].append(read[index2] + inser_diff * '-')
                else:
                    nalig['reads'][name].append(inser_diff * '-')
                    if just_once:
                        index2_delta += 1
                        just_once = False

        else:
            inser_diff = inser_2 - inser_1
            nalig['reference'].append(ref2_item)
            just_once = True
            for name, read in align1['reads'].items():
                if name not in nalig['reads']:
                    nalig['reads'][name] = []
                if inser_1 > 0:
                    nalig['reads'][name].append(read[index1] + inser_diff * '-')
                else:
                    nalig['reads'][name].append(inser_diff * '-')
                    if just_once:
                        index1_delta += 1
                        just_once = False
            add_collumn(align2, nalig, index)


    return nalig

def _make_multiple_alignment(alignments, reads=None):
    'It makes the multime alignments using ref to read sinmple alignments'
    orig_align = alignments
    snv_types_per_read = {}
    alignments = _prepare_alignments(alignments)
    if not alignments:
        return None
    alignment = alignments[0]
    assert len(alignment['reads'].values()[0]) == len(alignment['reference'])
    for index in range(1, len(alignments)):
        assert len(alignments[index]['reads'].values()[0]) == len(alignments[index]['reference'])
        alignment = _join_alignments(alignment, alignments[index],
                                     snv_types_per_read)
    try:
        assert len(alignment['reads'].values()[0]) == len(alignment['reference'])
    except AssertionError:
        print 'input alignments: ', orig_align
        print 'resul_alignment: ', alignment
        raise
    return alignment

def _get_alignment_section(pileup_read, start, end, reference_seq=None):
    'It gets a section of the alignment of the given read'
    # we don't get the last position because our politic is:
    #  [start, end[
    #  [start, stop]

    stop = end - 1

    aligned_read = pileup_read.alignment
    read_seq     = aligned_read.seq

    (ref_segments, read_segments, ref_limits, read_limits, segment_types,
    segment_lens) = _get_cigar_segments_from_aligned_read(aligned_read)

    #in which segment starts the section
    start_segment = _locate_segment(start, ref_segments, segment_lens, ref_limits)
    start_segment = start_segment[0] if start_segment is not None else None

    end_segment  = _locate_segment(end, ref_segments, segment_lens, ref_limits)
    end_segment  = end_segment[0]  if end_segment  is not None else None

    stop_segment = _locate_segment(stop, ref_segments, segment_lens, ref_limits)
    stop_segment = stop_segment[0] if stop_segment is not None else None

    #when does the read start and end in the reference coordinate system
    ref_start_limit, ref_end_limit = ref_limits
    read_start_in_ref = start if start >= ref_start_limit else ref_start_limit
    read_stop_in_ref  = stop  if stop  <= ref_end_limit   else ref_end_limit

    # we have to look if  the position of the end is in an insertion. For that
    # we can look the difference between the end_segment and
    # the origial_end_segment
    # If the difference is bigger than 1 is because there is an insertion
    if (end_segment is not None and stop_segment is not None and
        end_segment - stop_segment > 1):
        stop_segment += 1

    cum_ref_seq, cum_read_seq = '', ''
    #before alignment
    len_before_segment = ref_start_limit - start
    if reference_seq is None:
        ref_seq_before = UNKNOWN_NUCLEOTYDE * len_before_segment
    else:
        ref_seq_before = str(reference_seq.seq[ref_start_limit - len_before_segment:ref_start_limit])
    cum_ref_seq += ref_seq_before
    cum_read_seq += NO_NUCLEOTYDE * len_before_segment

    #in alignment
    ssegment = 0 if start_segment is None else start_segment
    esegment = len(ref_segments) - 1 if stop_segment is None else stop_segment
    for isegment in range(ssegment, esegment + 1):
        seg_type = segment_types[isegment]
        seg_len = segment_lens[isegment]
        ref_seg_start = ref_segments[isegment]
        # when the segment is an insert there is no start, we have to calculate
        if ref_seg_start is None:
            ref_seg_start = ref_segments[isegment + 1] - 1

        read_seg_start = read_segments[isegment]
        if isegment == ssegment:
            start_delta = read_start_in_ref - ref_seg_start
        else:
            start_delta = 0

        if isegment == esegment:
            end_delta = seg_len - read_stop_in_ref + ref_seg_start - 1
#            end_delta =  read_stop_in_ref - ref_seg_start + 1
        else:
            end_delta = 0

        if seg_type == INSERTION:
#            seg_ref_seq = DELETION_ALLELE * (seg_len - end_delta)
            seg_ref_seq = DELETION_ALLELE * seg_len
            read_start = read_seg_start + start_delta
            read_end = read_start + seg_len  #- end_delta
            seg_read_seq = read_seq[read_start: read_end]
        elif seg_type == DELETION or seg_type == SKIP:
            ref_start = ref_seg_start + start_delta
            ref_end = ref_start + seg_len - start_delta - end_delta

            if reference_seq is not None:
                seg_ref_seq = str(reference_seq.seq[ref_start: ref_end])
            else:
                seg_ref_seq = UNKNOWN_NUCLEOTYDE * (seg_len - start_delta - end_delta)
            read_start = None
            read_end = None
            seg_read_seq = DELETION_ALLELE * (seg_len - start_delta - end_delta)

        else:
            ref_start = ref_seg_start + start_delta
            ref_end = ref_start + seg_len - end_delta - start_delta
            if reference_seq is not None:
                seg_ref_seq = str(reference_seq.seq[ref_start: ref_end])
            else:
                seg_ref_seq = UNKNOWN_NUCLEOTYDE * (seg_len - end_delta - start_delta)
            read_start = read_seg_start + start_delta
            read_end = read_start + seg_len - end_delta - start_delta
            seg_read_seq = read_seq[read_start: read_end]

        cum_ref_seq += seg_ref_seq
        cum_read_seq += seg_read_seq

    #after alignment
    len_after_segment = stop - ref_limits[1]

    if reference_seq is None:
        ref_seq_after = UNKNOWN_NUCLEOTYDE * len_after_segment
    else:
        ref_seq_after = str(reference_seq.seq[ref_limits[1] + 1:ref_limits[1] + len_after_segment + 1])
    cum_ref_seq  += ref_seq_after
    cum_read_seq += NO_NUCLEOTYDE * len_after_segment

    assert len(cum_ref_seq) == len(cum_read_seq)
    return cum_ref_seq, cum_read_seq

def _get_alleles_from_read(ref_allele, ref_pos, pileup_read):
    '''It returns an allele from the read.

    It returns a list with the alleles in the given position.
    The returned allele can be an empty list if we're in a deletion.
    If the position holds an insertion it will return two alleles, the
    insertion and the nucleotide at that position.
    '''

    alleles = []
    aligned_read = pileup_read.alignment

    (ref_segments, read_segments, ref_limits, read_limits, segment_types,
    segment_lens) = _get_cigar_segments_from_aligned_read(aligned_read)
    located_segment = _locate_segment(ref_pos, ref_segments, segment_lens,
                                      ref_limits)
    if located_segment is None:
        return []
    else:
        segment_index, segment_pos = located_segment
    is_reverse = bool(aligned_read.is_reverse)

    if segment_types[segment_index] == MATCH:
        read_pos = _from_ref_to_read_pos(MATCH, ref_segments[segment_index],
                                         read_segments[segment_index], ref_pos)
        allele, qual = _get_raw_allele_from_read(aligned_read, read_pos)
        if allele != ref_allele:
            kind = SNP
        else:
            kind = INVARIANT

        alleles.append((allele, kind, qual, is_reverse))

        if segment_pos == IN_LAST_POS or segment_pos == IN_FIRST_AND_LAST:
            #Is there an insertion in the next position?
            next_read_pos = read_pos + 1
            allele, kind, qual = _get_insertion(segment_index, segment_types,
                                                next_read_pos, aligned_read,
                                                segment_lens)
            if kind is not None:
                alleles.append((allele, kind, qual, is_reverse))

    elif segment_types[segment_index] == DELETION:
        if (segment_pos == IN_FIRST_POS or segment_pos == IN_FIRST_AND_LAST or
            segment_pos == IN_LAST_POS):
            read_pos1, read_pos2 = _read_pos_around_del(ref_segments,
                                                       read_segments,
                                                       segment_types,
                                                       segment_index,
                                                       ref_pos)

        if segment_pos == IN_FIRST_POS or segment_pos == IN_FIRST_AND_LAST:
            indel_length = segment_lens[segment_index]
            allele = DELETION_ALLELE * (indel_length)
            #in the deletion case the quality is the lowest of the
            #bases that embrace the deletion
            if aligned_read.qual:
                qual0 = aligned_read.qual[read_pos1]
                qual0 = _quality_to_phred(qual0)
                qual1 = aligned_read.qual[read_pos2]
                qual1 = _quality_to_phred(qual1)
                qual = min((qual0, qual1))
            else:
                qual = None
            kind = DELETION
            alleles.append((allele, kind, qual, is_reverse))
        if segment_pos ==IN_FIRST_AND_LAST or segment_pos == IN_LAST_POS:
            #Is there an insertion in the next position?
            allele, kind, qual = _get_insertion(segment_index, segment_types,
                                               read_pos1, aligned_read,
                                               segment_lens)
            if kind is not None:
                alleles.append((allele, kind, qual, is_reverse))

    elif segment_types[segment_index] == SKIP:
        #Is there an insertion in the next position?
        read_pos1, read_pos2 = _read_pos_around_del(ref_segments,
                                                       read_segments,
                                                       segment_types,
                                                       segment_index,
                                                       ref_pos)
        allele, kind, qual = _get_insertion(segment_index, segment_types,
                                           read_pos1, aligned_read,
                                           segment_lens)
        if kind is not None:
            alleles.append((allele, kind, qual, is_reverse))

    elif segment_types[segment_index] == INSERTION:
        pass    #if we're in an insertion, it is returned in the last position
                #of the previous match segment
    return alleles, read_limits

def _quality_to_phred(quality):
    'It transforms a qual chrs into a phred quality'
    if quality is None:
        return None
    elif len(quality) == 1:
        phred_qual = ord(quality) - 33
    else:
        phred_quals = [ord(qual) - 33 for qual in quality]
        phred_qual = sum(phred_quals) / len(phred_quals)
    if phred_qual == 93:  #the character used for unknown qualities
        phred_qual = None
    return phred_qual

def _add_allele(alleles, allele, kind, read_name, read_group, is_reverse, qual,
                mapping_quality, readgroup_info, pileup_read):
    'It adds one allele to the alleles dict'
    key = (allele, kind)
    if key not in alleles:
        alleles[key] = {'read_groups':[], 'orientations':[],
                        'qualities':[], 'mapping_qualities':[], 'reads':[]}
    allele_info = alleles[key]
    allele_info['read_groups'].append(read_group)
    allele_info['orientations'].append(not(is_reverse))
    allele_info['qualities'].append(qual)
    allele_info['mapping_qualities'].append(mapping_quality)
    allele_info['reads'].append(pileup_read)

def _normalize_read_edge_conf(read_edge_conf):
    'It returns a dict with all valid keys'
    platforms = ('454', 'sanger', 'illumina')
    if read_edge_conf is None:
        read_edge_conf = {}
    for platform in platforms:
        if platform not in read_edge_conf:
            read_edge_conf[platform] = (None, None)
    return read_edge_conf

def _update_pileup_reads(snv, reads):
    'It gets the pileup read for each snv_allele'
    snv_name = snv['ref_name'] + '_' + str(snv['ref_position'])
    if snv_name not in reads:
        reads[snv_name] = {}
    for allele, allele_info in snv['alleles'].items():
        reads[snv_name][allele] = allele_info['reads']


def check_read_length(read, position):
    'If checks if the given pileup read fills all the snv extension'
    aligned_read = read.alignment
    read_in_ref = _get_cigar_segments_from_aligned_read(aligned_read)[2]
    if read_in_ref[0] < position[0] and read_in_ref[1] > position[1]:
        return True
    else:
        return False

def _check_enough_reads(reads, position):
    '''It checks that the reads in the snv are of the length of the snv.
    Otherwise they will be removed'''

    for snv_name, alleles in reads.items():
        for allele, pileup_reads in alleles.items():
            new_reads = []
            for read in pileup_reads:
                if check_read_length(read, position):
                    new_reads.append(read)
            if not new_reads:
                del reads[snv_name][allele]
            else:
                reads[snv_name][allele] = new_reads
        alleles = reads[snv_name].keys()
        if not alleles or (len(alleles) == 1 and alleles[0][1] == INVARIANT):
            del reads[snv_name]
            return False
    return True

def _increase_snv_length(snv, snv_block, reads):
    '''It checks if we have increase the snv or not. It will check if the reads
    have the proper length and that the increase will produce by deletion
    elongation'''

    if snv_block['end'] <= _get_snv_start_vcf4(snv):
        return False
    posible_end = _get_snv_end_position(snv)
    snv_pos = (snv_block['start'], posible_end)
    return _check_enough_reads(reads, snv_pos)

def _get_snv_start_vcf4(snv):
    '''It corrects the start of the snv to the new vcf4 system

           01234567      snp_caller                 vcf_format
      ref  atctgtag
    read1  atccgtag       snp in pos 3             snp in pos 3
    read2  atctgcctag     inser in pos 5           inser in pos 5
    read3  atgcctag       del in pos 2             del in pos 1

    This is the result and
    '''

    if (_is_snv_of_kind(snv, DELETION) or _is_snv_of_kind(snv, INDEL) or
        _is_snv_of_kind(snv, COMPLEX)):
        start = snv['ref_position'] - 1
    else:
        start = snv['ref_position']
    return start

def _get_snv_end_position(snv):
    'it returns the snv position plus the length of the allele'

    end = snv['ref_position']
    if _is_snv_of_kind(snv, INSERTION) or _is_snv_of_kind(snv, SNP):
        end += 1
    else:
        max_length = 0
        for allele in snv['alleles'].keys():
            allele, kind = allele
            allele_length = len(allele)
            if kind == DELETION and allele_length > max_length:
                max_length = allele_length
        end += max_length
    return end

def _make_snv_blocks(snvs):
    '''It joins snvs that should be just one snv. e.g. a deletion that match
    with another deletion in the same position
    '''
    snv_block = {}
    reads     = {}
    snvs1, snvs2 = itertools.tee(snvs, 2)
    try:
        snvs2.next() # pass the first one
    except StopIteration:
        pass
    for snv in snvs1:
#        if  snv['ref_position'] == 724056:
#            print snv['alleles'].keys()
#            pass
        try:
            next_snv = snvs2.next()
        except StopIteration:
            next_snv = None

        start = _get_snv_start_vcf4(snv)
        end   = _get_snv_end_position(snv)
        if not snv_block:
            snv_block = {'start':start, 'end':end, 'snvs':[snv]}
        else:
            if _increase_snv_length(snv, snv_block, reads):
                if snv_block['end'] < end:
                    snv_block['end'] = end
                snv_block['snvs'].append(snv)
            else:
                if (next_snv is not None and
                    _get_snv_start_vcf4(snv) == start + 1 and
                    _is_snv_of_kind(next_snv, DELETION)):
                    if _increase_snv_length(next_snv, snv_block, reads):
                        new_snv_end = _get_snv_end_position(next_snv)
                        if snv_block['end'] < new_snv_end:
                            snv_block['end'] = new_snv_end
                        snv_block['snvs'].append(next_snv)

                    # don't use the next snp as we have checked in this step
                    snvs1.next()
                    try:
                        snvs2.next()
                    except StopIteration:
                        pass
                else:
                    yield snv_block
                    reads = {}
                    snv_block = {'start':start, 'end':end, 'snvs':[snv]}
        _update_pileup_reads(snv, reads)
    else:
        if snv_block:
            yield snv_block

def get_insertions_in_position(reads, position):
    'it returns the insertions that are found inside de positions'
    insertions = set()
    for segments in reads.values():
        (ref_segments, read_segments, ref_limits, read_limits, segment_types,
        segment_lens) = segments['segments']
        for index, segment_type in enumerate(segment_types):
            if segment_type == 1:
                insert_start_in_ref = ref_segments[index - 1]
                insert_end_in_ref = insert_start_in_ref +  segment_lens[index]
                if ((insert_start_in_ref > position[0] and insert_start_in_ref< position[1]) or
                    (insert_end_in_ref > position[0] and insert_end_in_ref< position[1])):
                    insertions.add((insert_start_in_ref, insert_end_in_ref))

    return list(insertions)

def _is_snv_of_kind(snv, kind):
    'True if it is a insertion false if not'
    allele_types = [allele[1] for allele in snv['alleles'].keys()]
    snv_kind = _calculate_snv_kinds(allele_types)
    if snv_kind == kind:
        return True
    else:
        return False

def _calculate_snv_kinds(kinds):
    'It returns the snv kind for the given feature'
    if len(kinds) == 1:
        return kinds[0]
    kind = kinds[0]
    for index in range(1, len(kinds)):
        kind = _calculate_kind(kind, kinds[index])
    return kind

def _calculate_allele_kind(ref_allele, allele):
    'It calculates the allele kind'
    if ref_allele == allele:
        return INVARIANT
    kinds = []
    for ref_n, allele_n in zip(ref_allele, allele):
        if ref_n == allele_n or ref_n == '' or allele_n == '':
            kind = INVARIANT
        elif '-' in ref_n :
            kind = INSERTION
        elif '-' in allele_n:
            kind = DELETION
        else:
            kind = SNP
        kinds.append(kind)
    return _calculate_snv_kinds(kinds)


def _remove_allele_from_alignments(alignments, min_num_reads_for_allele):
    ''' It removes alignments/alleles taking into account the times it appears
    min_num_reads_for_allele'''
    allele_count = {}
    for ref, allele in alignments.values():
        if allele not in allele_count:
            allele_count[allele] = 0
        allele_count[allele] +=1

    for read_name, alignment in alignments.items():
        read_allele = alignment[1]
        if allele_count[read_allele] < min_num_reads_for_allele:
            del alignments[read_name]


def _join_snvs(snv_block, min_num_alleles, min_num_reads_for_allele,
               min_quality, reference_seq=None):
    'It joins the snvs that should be together'
    snvs = snv_block['snvs']
    if len(snvs) == 1 and _is_snv_of_kind(snvs[0], SNP):
        yield snvs[0]
    else:
        start = snv_block['start']
        end   = snv_block['end']
        new_snv = {'ref_name':snvs[0]['ref_name'],
                   'ref_position':start ,
                   'read_groups':snvs[0]['read_groups'],
                   'alleles':{}}

        # collect all the reads and its alignment
        reads = {}
        alignments = {}
        allele_kinds = {}
        for snv in snvs:
            for allele, allele_info in snv['alleles'].items():
                allele_kind = allele[1]
                for index in range(len(allele_info['reads'])):
                    name = allele_info['reads'][index].alignment.qname
                    if name in reads:
                        continue
                    reads[name] = {}
                    read        = allele_info['reads'][index]
                    read_group  = allele_info['read_groups'][index]
                    orientation = allele_info['orientations'][index]
                    quality     = allele_info['qualities'][index]
                    mapping_quality =  allele_info['mapping_qualities'][index]

                    reads[name]['read'] = read
                    reads[name]['read_group'] = read_group
                    reads[name]['orientation'] = orientation
                    reads[name]['quality']= quality
                    reads[name]['mapping_quality']= mapping_quality

                for index, read in enumerate(allele_info['reads']):
                    read_name = read.alignment.qname
                    if read_name not in allele_kinds:
                        allele_kinds[read_name] = []
                    allele_kinds[read_name].append(allele_kind)
                    if read_name not in alignments:
                        alignment = _get_alignment_section(read, start, end,
                                                    reference_seq=reference_seq)
#                        print read_name, alignment, start, end
                        alignments[read_name] = alignment
        _remove_allele_from_alignments(alignments, min_num_reads_for_allele)
        try:
            malignment = _make_multiple_alignment(alignments, reads)
#            alis = []
#            for ali in malignment['reads'].values():
#                if ali not in alis:
#                    alis.append(ali)
#            print 'alis', alis
        except AssertionError:
            print 'ref_name: ', snvs[0]['ref_name']
            print 'position: ', start
            raise

        if malignment is not None:
            ref_allele = malignment['reference']
            for read, allele in malignment['reads'].items():
                kind = _calculate_allele_kind(ref_allele, allele)

                allele = ''.join(allele)
                allele = (allele, kind)
                if allele not in new_snv['alleles']:
                    new_snv['alleles'][allele] = {'read_groups':[],
                                                  'reads':[],
                                                  'orientations':[],
                                                  'qualities':[],
                                                  'mapping_qualities':[]}
                new_snv['alleles'][allele]['read_groups'].append(reads[read]['read_group'])
                new_snv['alleles'][allele]['reads'].append(reads[read]['read'])
                new_snv['alleles'][allele]['orientations'].append(reads[read]['orientation'])
                new_snv['alleles'][allele]['qualities'].append(reads[read]['quality'])
                new_snv['alleles'][allele]['mapping_qualities'].append(reads[read]['mapping_quality'])
        else:
            new_snv = None

        if new_snv is not None:
            # calculate snv quality
            for allele_info in new_snv['alleles'].values():
                allele_info['quality'] = _calculate_allele_quality(allele_info)
            #reference allele

            new_snv['reference_allele'] = ''.join(ref_allele)

            alleles = new_snv['alleles']

            #remove bad quality alleles
            _remove_bad_quality_alleles(alleles, min_quality)

            # min_num_reads_for_allele
            _remove_alleles_by_read_number(alleles, min_num_reads_for_allele)

            #if there are a min_num number of alleles requested and there are more
            #alleles than that
            #OR
            #there is some allele different than invariant
            #a variation is yield
            if (len(alleles) > min_num_alleles or
                (min_num_alleles == 1 and alleles.keys()[0][1] != INVARIANT) or
                (min_num_alleles > 1 and len(alleles) >= min_num_alleles)):
                yield new_snv



def _snvs_in_bam(bam, reference, min_quality,
                             default_sanger_quality, min_mapq, min_num_alleles,
                             max_maf, min_num_reads_for_allele,
                             read_edge_conf=None, default_bam_platform=None):
    'It returns the snvs in a bam for the given reference'
    snvs = _snvs_in_bam_by_position(bam, reference, min_quality,
                                default_sanger_quality, min_mapq,
                                min_num_alleles,max_maf,
                                min_num_reads_for_allele, read_edge_conf,
                                default_bam_platform)
    for snv_block in _make_snv_blocks(snvs):
        for snv in _join_snvs(snv_block, min_num_alleles,
                              min_num_reads_for_allele, min_quality,
                              reference_seq=reference):
            yield snv

def _snvs_in_bam_by_position(bam, reference, min_quality,
                             default_sanger_quality, min_mapq, min_num_alleles,
                             max_maf, min_num_reads_for_allele,
                             read_edge_conf, default_bam_platform):
    '''It yields the snv information for every snv in the given reference,
    for each position'''

    min_num_alleles = int(min_num_alleles)

    read_groups_info = get_read_group_info(bam)
    if not read_groups_info:
        if default_bam_platform is None:
            msg = 'Platform is not present either in header or in '
            msg += 'configuration'
            raise ValueError(msg)
        read_groups_info = {UNKNOWN_RG:{'PL':default_bam_platform}}

    reference_id = get_seq_name(reference)
    reference_seq = reference.seq
    reference_len = len(reference_seq)
    #we can clean the cache of segments because we're in a new molecule
    global SEGMENTS_CACHE
    SEGMENTS_CACHE = {}
    for column in bam.pileup(reference=reference_id):
        alleles = {}
        ref_pos = column.pos
        if ref_pos >= reference_len:
            continue
        ref_id = bam.getrname(column.tid)
        ref_allele = reference_seq[ref_pos].upper()
        for pileup_read in column.pileups:
            #for each read in the column we add its allele to the alleles dict
            aligned_read = pileup_read.alignment

            read_mapping_qual = aligned_read.mapq
            #We ignore the reads that are likely to be missaligned
            if read_mapping_qual < min_mapq:
                continue

            try:
                read_group = aligned_read.opt('RG')
            except KeyError:
                read_group = UNKNOWN_RG

            read_name = aligned_read.qname
            if read_groups_info and read_group in read_groups_info:
                platform = read_groups_info[read_group]['PL']
            else:
                platform = default_bam_platform

            read_pos = pileup_read.qpos


            alleles_here, read_limits = _get_alleles_from_read(ref_allele,
                                                               ref_pos,
                                                               pileup_read)
            #if read_name == '964_643_534_F3':
            #    print alleles_here, read_pos, ref_pos

            if read_edge_conf and platform in read_edge_conf:
                edge_left, edge_right = read_edge_conf[platform]

                #if we're in the edge region to be ignored we continue to
                #the next read, because there's no allele to add for this one.

                if (edge_left is not None and
                    read_limits[0] + edge_left > read_pos):
                    continue
                if (edge_right is not None and
                    read_pos > read_limits[1] - edge_right):
                    continue

            for allele in alleles_here:
                allele, kind, qual, is_reverse = allele
                _add_allele(alleles, allele, kind, read_name, read_group,
                    is_reverse, qual, read_mapping_qual,
                    read_groups_info, pileup_read)

        #remove N
        _remove_alleles_n(alleles)

        #add default sanger qualities to the sanger reads with no quality
        _add_default_sanger_quality(alleles, default_sanger_quality,
                                    read_groups_info)
        if ref_pos == 724058:
            print alleles.keys()

        #remove bad quality alleles
        _remove_bad_quality_alleles(alleles, min_quality)

        #check maf
        if not check_maf_ok(alleles, max_maf):
            continue

        # min_num_reads_for_allele
        _remove_alleles_by_read_number(alleles, min_num_reads_for_allele)

        #if there are a min_num number of alleles requested and there are more
        #alleles than that
        #OR
        #there is some allele different than invariant
        #a variation is yield
        if not alleles:
            continue
        if (len(alleles) > min_num_alleles or
            (min_num_alleles == 1 and alleles.keys()[0][1] != INVARIANT) or
            (min_num_alleles > 1 and len(alleles) >= min_num_alleles)):
            yield {'ref_name':ref_id,
                   'ref_position':ref_pos,
                   'reference_allele':ref_allele,
                   'alleles':alleles,
                   'read_groups':read_groups_info}

def _remove_alleles_by_read_number(alleles, min_num_reads_for_allele):
    'It remove alleles with less reads than the given value'
    alleles_to_remove = []
    for allele_name, allele_info in  alleles.items():
        if len(allele_info['read_groups']) < min_num_reads_for_allele:
            alleles_to_remove.append(allele_name)

    if alleles_to_remove:
        for allele_to_remove in alleles_to_remove:
            del(alleles[allele_to_remove])

def _add_default_sanger_quality(alleles, default_sanger_quality,
                                read_groups_info):
    'It adds default sanger qualities to the sanger reads with no quality'

    for allele_info in alleles.values():
        for index, (qual, rg) in enumerate(zip(allele_info['qualities'],
                                             allele_info['read_groups'])):
            try:
                if qual is None and read_groups_info[rg]['PL'] == 'sanger':
                    allele_info['qualities'][index] = default_sanger_quality
            except KeyError:
                if 'PL' not in read_groups_info[rg]:
                    msg = 'The bam file has no platforms for the read groups'
                    raise KeyError(msg)
                else:
                    raise

def _remove_alleles_n(alleles):
    'It deletes the aleles that are N'
    for allele in alleles:
        if allele[0] in N_ALLLELES:
            del alleles[allele]

def _remove_bad_quality_alleles(alleles, min_quality):
    'It adds the quality to the alleles dict and it removes the bad alleles'

    orientations_independent = False
    if orientations_independent:
        qual_calculator = _calculate_allele_quality_oriented
    else:
        qual_calculator = _calculate_allele_quality
    for allele, allele_info in alleles.items():
        qual = qual_calculator(allele_info)
        allele_info['quality'] = qual
        if qual < min_quality:
            del alleles[allele]

def _calculate_allele_quality(allele_info):
    'It returns the quality for the given allele'

    #we sort all qualities
    quals = allele_info['qualities'][:]
    #slow alternative
    #quals.sort(lambda x, y: int(y - x))
    #fast alternative
    qual_set = set(quals)
    for index in range(3):
        if not qual_set:
            break
        quals[index] = max(qual_set)
        qual_set.remove(quals[index])
    total_qual = 0
    if quals:
        total_qual += quals[0]
        if len(quals) > 1:
            total_qual += quals[1] / 4.0
            if len(quals) > 2:
                total_qual += quals[2] / 4.0
    return total_qual

def _calculate_allele_quality_oriented(allele_info):
    '''It returns the quality for the given allele
    It assumes that reads with different orientations are independent'''
    #we gather all qualities for independent groups
    quals = defaultdict(list)
    for qual, orientation in zip(allele_info['qualities'],
                                 allele_info['orientations']):
        quals[orientation].append(qual)

    #we sort all qualities
    for independent_quals in quals.values():
        independent_quals.sort(lambda x, y: int(y - x))

    total_qual = 0
    for independent_quals in quals.values():
        if independent_quals:
            total_qual += independent_quals[0]
            if len(independent_quals) > 1:
                total_qual += independent_quals[1] / 4.0
                if len(independent_quals) > 2:
                    total_qual += independent_quals[2] / 4.0
    return total_qual

def _root_mean_square(numbers):
    'It returns the root mean square for the given numbers'
    power2 = lambda x: math.pow(x, 2)
    return math.sqrt(sum(map(power2, numbers)) / len(numbers))

def _summarize_snv(snv):
    'It returns an snv with an smaller memory footprint'
    used_read_groups = set()
    for allele_info in snv['alleles'].values():
        #the read_groups list to a count dict
        rg_count = {}
        for read_group in allele_info['read_groups']:
            if read_group not in rg_count:
                rg_count[read_group] = 0
            rg_count[read_group] += 1
            used_read_groups.add(read_group)
        allele_info['read_groups'] = rg_count

    #we calculate a couple of parameters that summarize the quality
    for kind in ('mapping_qualities', 'qualities'):
        quals = []
        for allele_info in snv['alleles'].values():
            quals.extend(allele_info[kind])
        if kind == 'mapping_qualities':
            kind = 'mapping_quality'
        if kind == 'qualities':
            kind = 'quality'
        snv[kind] = _root_mean_square(quals) if quals else None

    for allele_info in snv['alleles'].values():
        #we remove some extra quality info
        del allele_info['mapping_qualities']
        del allele_info['qualities']
        del allele_info['orientations']
        del allele_info['reads']

    #we remove from the read_groups the ones not used in this snv
    new_read_groups = {}
    for read_group, info in snv['read_groups'].items():
        if read_group in used_read_groups:
            new_read_groups[read_group] = info

    snv['read_groups'] = new_read_groups
    return snv

def create_snv_annotator(bam_fhand, min_quality=45, default_sanger_quality=25,
                         min_mapq=15, min_num_alleles=1, max_maf=None,
                         read_edge_conf=None, default_bam_platform=None,
                         min_num_reads_for_allele=None, ploidy=2):
    'It creates an annotator capable of annotating the snvs in a SeqRecord'

    #the bam should have an index, does the index exists?
    bam_fhand = get_fhand(bam_fhand)
    create_bam_index(bam_fpath=bam_fhand.name)
    read_edge_conf = _normalize_read_edge_conf(read_edge_conf)

    bam = pysam.Samfile(bam_fhand.name, 'rb')
    # default min num_reads per allele and ploidy
    if min_num_reads_for_allele is None:
        min_num_reads_for_allele = DEFAUL_MIN_NUM_READS_PER_ALLELE
    if ploidy is None:
        ploidy = DEFAULT_PLOIDY

    def annotate_snps(sequence):
        'It annotates the snvs found in the sequence'
        for snv in _snvs_in_bam(bam, reference=sequence,
                                min_quality=min_quality,
                                default_sanger_quality=default_sanger_quality,
                                min_mapq=min_mapq,
                                min_num_alleles=min_num_alleles,
                                max_maf=max_maf,
                                read_edge_conf=read_edge_conf,
                                default_bam_platform=default_bam_platform,
                             min_num_reads_for_allele=min_num_reads_for_allele):
            snv = _summarize_snv(snv)
            location = snv['ref_position']
            type_ = 'snv'

            qualifiers = {'alleles':snv['alleles'],
                          'reference_allele':snv['reference_allele'],
                          'read_groups':snv['read_groups'],
                          'mapping_quality': snv['mapping_quality'],
                          'quality': snv['quality']}
            snv_feat = SeqFeature(location=FeatureLocation(location, location),
                              type=type_,
                              qualifiers=qualifiers)

            annotate_pic(snv_feat)
            annotate_heterozygosity(snv_feat, ploidy=ploidy)

            sequence.features.append(snv_feat)
        return sequence
    return annotate_snps

def calculate_snv_kind(feature, detailed=False):
    'It returns the snv kind for the given feature'
    snv_kind = INVARIANT
    alleles = feature.qualifiers['alleles']
    for allele in alleles.keys():
        allele_kind = allele[1]
        snv_kind = _calculate_kind(allele_kind, snv_kind)

    if snv_kind == SNP and detailed:
        snv_kind = _guess_snp_kind(alleles)

    return snv_kind

def _al_type(allele):
    'I guesses the type of the allele'
    allele = allele.upper()
    if allele in ('A', 'G'):
        return 'purine'
    elif allele in ('T', 'C'):
        return 'pirimidine'
    return UNKNOWN

def _guess_snp_kind(alleles):
    'It guesses the type of the snp'
    alleles = alleles.keys()
    # if we take into account the reference to decide if there is a variation
    if len(alleles) < 2:
        return UNKNOWN
    al0 = _al_type(alleles[0][0])
    al1 = _al_type(alleles[1][0])
    if al0 == UNKNOWN or al1 == UNKNOWN:
        snv_kind = UNKNOWN
    elif al0 == al1:
        snv_kind = TRANSITION
    else:
        snv_kind = TRANSVERSION
    return snv_kind

def _calculate_kind(kind1, kind2):
    'It calculates the result of the union of two kinds'
    if kind1 == kind2:
        return kind1
    else:
        if kind1 is INVARIANT:
            return kind2
        elif kind2 is INVARIANT:
            return kind1
        elif kind1 in [SNP, COMPLEX] or kind2 in [SNP, COMPLEX]:
            return COMPLEX
        else:
            return INDEL

def _cmp_by_read_num(allele1, allele2):
    'cmp by the number of reads for each allele'
    return len(allele2['read_names']) - len(allele1['read_names'])

def sorted_alleles(feature):
    'It returns the alleles sorted by number of reads'
    #from dict to list
    alleles = feature.qualifiers['alleles']
    alleles_list = []
    for allele, allele_info in alleles.items():
        allele_info = copy(allele_info)
        allele_info['seq'] = allele[0]
        allele_info['kind'] = allele[1]
        alleles_list.append(allele_info)
    return sorted(alleles_list, _cmp_by_read_num)

def snvs_in_window(snv, snvs, window, snv_type=None, maf=None):
    'it gets all the snvs in a window taking a snv as reference'

    num_of_snvs = 0
    snv_location = int(str(snv.location.start))
    left_margin = snv_location - (window / 2)
    rigth_margin = snv_location + (window / 2)
    for snv in snvs:
        current_location = int(str(snv.location.start))
        if current_location == snv_location:
            continue
        if current_location >= left_margin and current_location <= rigth_margin:

            if snv_type is None and maf is None:
                num_of_snvs += 1
            elif snv_type is None and maf is not None:
                snv_maf = calculate_maf_frequency(snv)
                if snv_maf < maf:
                    num_of_snvs += 1
            elif snv_type is not None and maf is None:
                type_ = calculate_snv_kind(snv)
                if ((snv_type == type_) or
                    (snv_type == INDEL and type_ in(INSERTION, DELETION))):
                    num_of_snvs += 1
            else:
                type_ = calculate_snv_kind(snv)
                snv_maf = calculate_maf_frequency(snv)
                if (((snv_type == type_) or
                    (snv_type == INDEL and type_ in(INSERTION, DELETION))) and
                    (snv_maf < maf)):
                    num_of_snvs += 1
    return num_of_snvs

def _get_group(read_group, group_kind, read_groups):
    'It returns the group (lb, rg, sm) for the given rg and group_kind'
    if group_kind:
        if group_kind == 'read_groups':
            return read_group
        else:
            group_kind = group_kind.lower()
            if group_kind in ('lb', 'library', 'libraries'):
                group_kind = 'LB'
            elif group_kind in ('sm', 'sample', 'samples'):
                group_kind = 'SM'
            elif group_kind in ('pl', 'platform', 'platforms'):
                group_kind = 'PL'
            return read_groups[read_group][group_kind]

def check_maf_ok(alleles, max_maf):
    'It checks that the major allele freq is less than maximun limit'
    if max_maf is None:
        return True
    maf = _calculate_maf_frequency_for_alleles(alleles,alleles_is_dict=False)
    if maf > max_maf:
        return False
    else:
        return True

def _allele_count(allele, alleles, read_groups=None,
                  groups=None, group_kind=None, alleles_is_dict=True):
    'It returns the number of reads for the given allele'
    counts = []
    if not alleles_is_dict:
        return len(alleles[allele]['read_groups'])
    for read_group, count in alleles[allele]['read_groups'].items():
        #do we have to count this read_group?
        group = _get_group(read_group, group_kind, read_groups)
        if not groups or groups and group in groups:
            counts.append(count)
    return sum(counts)

def calculate_maf_frequency(feature, groups=None, group_kind=None):
    'It returns the most frequent allele frequency'
    alleles     = feature.qualifiers['alleles']
    read_groups = feature.qualifiers['read_groups']

    return _calculate_maf_frequency_for_alleles(alleles, groups=groups,
                                                group_kind=group_kind,
                                                read_groups=read_groups)

def _calculate_maf_frequency_for_alleles(alleles, groups=None, group_kind=None,
                                        read_groups=None, alleles_is_dict=True):
    'It returns the most frequent allele frequency'
    major_number_reads = None
    total_number_reads = 0
    for allele in alleles:
        number_reads = _allele_count(allele, alleles, read_groups, groups,
                                     group_kind, alleles_is_dict)
        if major_number_reads is None or major_number_reads < number_reads:
            major_number_reads = number_reads
        total_number_reads += number_reads
    if not total_number_reads:
        return None
    return major_number_reads / total_number_reads

def calculate_snv_variability(sequence):
    'It returns the number of snv for every 100 pb'
    n_snvs = sum(1 for snv in sequence.get_features(kind='snv'))
    return n_snvs / len(sequence)

def calculate_cap_enzymes(feature, sequence, all_enzymes=False):
    '''Given an snv feature and a sequence it returns the list of restriction
    enzymes that distinguish between their alleles.'''

    if 'cap_enzymes' in feature.qualifiers:
        return feature.qualifiers['cap_enzymes']

    #which alleles do we have?
    alleles = set()
    for allele in feature.qualifiers['alleles'].keys():
        alleles.add(repr((allele[0], allele[1])))
    #for every pair of different alleles we have to look for differences in
    #their restriction maps
    enzymes = set()
    alleles = list(alleles)
    reference = sequence
    location = int(str(feature.location.start))
    for i_index in range(len(alleles)):
        for j_index in range(i_index, len(alleles)):
            if i_index == j_index:
                continue
            allelei = eval(alleles[i_index])
            allelei = {'allele':allelei[0], 'kind':allelei[1]}
            allelej = eval(alleles[j_index])
            allelej = {'allele':allelej[0], 'kind':allelej[1]}
            i_j_enzymes = _cap_enzymes_between_alleles(allelei, allelej,
                                                       reference, location,
                                                       all_enzymes)
            enzymes = enzymes.union(i_j_enzymes)

    enzymes = [str(enzyme) for enzyme in enzymes]
    feature.qualifiers['cap_enzymes'] = enzymes
    return enzymes

def _cap_enzymes_between_alleles(allele1, allele2, reference, location,
                                 all_enzymes=False):
    '''It looks in the enzymes that differenciate the given alleles.

    It returns a set.
    '''
    kind1 = allele1['kind']
    kind2 = allele2['kind']
    allele1 = allele1['allele']
    allele2 = allele2['allele']

    #we have to build the two sequences
    if all_enzymes:
        restriction_batch = CommOnly
    else:
        restriction_batch = RestrictionBatch(COMMON_ENZYMES)

    seq1 = create_alleles('seq1', allele1, kind1, reference, location)
    seq2 = create_alleles('seq2', allele2, kind2, reference, location)

    anal1 = Analysis(restriction_batch, seq1, linear=True)
    enzymes1 = set(anal1.with_sites().keys())
    anal1 = Analysis(restriction_batch, seq2, linear=True)
    enzymes2 = set(anal1.with_sites().keys())

    enzymes = set(enzymes1).symmetric_difference(set(enzymes2))

    return enzymes

def create_alleles(name, allele, kind, ref, loc):
    'The returns the sequence for the given allele'
    sseq = ref.seq
    if kind == INVARIANT:
        seq = sseq
    elif kind == SNP:
        seq = sseq[0:loc] + allele + sseq[loc + 1:]
    elif kind == DELETION:
        seq = sseq[0:loc + 1] + sseq[loc + len(allele) + 1:]
    elif kind == INSERTION:
        seq = sseq[0:loc] + allele + sseq[loc:]
    return seq

def variable_in_groupping(snv, group_kind, groups, in_union=True,
                           in_all_groups=True, reference_free=True, maf=None,
                           min_num_reads=None, min_reads_per_allele=None):
    'It looks if the given snv is variable for the given groups'

    alleles = _get_alleles_for_group(snv.qualifiers['alleles'],
                                     groups, group_kind,
                                     snv.qualifiers['read_groups'],
                                     min_reads_per_allele=min_reads_per_allele)
    if not alleles:
        return False

    if in_union:
        alleles = _aggregate_alleles(alleles)

    variable_per_read_group = []
    for alleles_in_rg in alleles.values():
        variable = _is_rg_variable(alleles_in_rg, reference_free=reference_free,
                                       maf=maf, min_num_reads=min_num_reads)
        variable_per_read_group.append(variable)

    if in_all_groups:
        return all(variable_per_read_group)
    else:
        return any(variable_per_read_group)

def invariant_in_groupping(snv, group_kind, groups, in_union=True,
                           in_all_groups=True, reference_free=True, maf=None,
                           min_num_reads=None):
    'it check if the given snv is invariant form the given groups'
    alleles = _get_alleles_for_group(snv.qualifiers['alleles'],
                                     groups, group_kind,
                                     snv.qualifiers['read_groups'])
    if not alleles and reference_free:
        return False
    elif not alleles and not reference_free:
        return True

    if in_union:
        alleles = _aggregate_alleles(alleles)

    invariable_per_read_group = []
    for alleles_in_rg in alleles.values():
        invariable = not _is_rg_variable(alleles_in_rg,
                                         reference_free=reference_free, maf=maf,
                                         min_num_reads=min_num_reads)
        invariable_per_read_group.append(invariable)

    if in_all_groups:
        return all(invariable_per_read_group)
    else:
        return any(invariable_per_read_group)

def _is_rg_variable(alleles, reference_free=True, maf=None, min_num_reads=None):
    'It checks if the allele is variable'
    allele_count = len(alleles)

    ref_in_alleles = False
    for allele in alleles.keys():
        if allele[1] == INVARIANT:
            ref_in_alleles = True

    if allele_count == 1:
        if reference_free:
            return False
        elif not reference_free and ref_in_alleles:
            return False

    maf_allele, num_reads = calc_maf_and_num_reads(alleles)
    if ((maf and maf < maf_allele) or
        (min_num_reads and min_num_reads > num_reads)):
        return False

    return True

def calc_maf_and_num_reads(alleles):
    'It calculates the maf and the number of reads in that position'
    values = alleles.values()
    major = max(values)
    num_reads  = sum(values)
    maf_allele = major/float(num_reads)
    return  maf_allele, num_reads

def _aggregate_alleles(alleles):
    'It joins all alleles for the read groups into one'
    aggregate = {}
    for alleles in alleles.values():
        for allele, allele_count in alleles.items():
            if allele not in aggregate:
                aggregate[allele] = 0
            aggregate[allele] += allele_count
    return {None: aggregate}

def _get_alleles_for_group(alleles, groups, group_kind='read_groups',
                           read_groups=None, min_reads_per_allele=None):
    '''It gets the alleles from the given items of type:key, separated by items.
    For example, if you give key rg and items rg1, rg2, it will return
    alleles separated in rg1 and rg2 '''
    alleles_for_groups = {}
    for allele, alleles_info in alleles.items():
        #print allele, alleles_info
        for read_group in alleles_info['read_groups']:
            if (min_reads_per_allele and
                alleles_info['read_groups'][read_group] < min_reads_per_allele):
                continue
            group = _get_group(read_group, group_kind, read_groups)
            if group not in groups:
                continue
            if not group in alleles_for_groups:
                alleles_for_groups[group] = {}
#            if allele not in alleles_for_groups[group]:
#                alleles_for_groups[group][allele] = 0
            alleles_for_groups[group][allele] = alleles_info['read_groups'][read_group]
    return alleles_for_groups

def calculate_pic(snv, read_groups=None, group_kind=None, groups=None):
    '''It calculates the uniformly minimum variance unbiased (UMVU) estimator
    of PIC of a locus, given a list with the number of times that each allele
    has been read.

    PIC(UMVU) = 1 - summatory((Xi(Xi-1))/(n(n-1))) -
                summatory((Xi(Xi-1)Xj(Xj-1))/(n(n-1)(n-2)(n-3))

    Xi = number of times that allele i-th has been read
    Xj = number of times that allele j-th has been read
    n = total number of reads

    Formula taken from "On Estimating the Heterozygosity and Polymorphism
    Information Content Value" by Shete S., Tiwari H. and Elston R.C.
    Theoretical Population Biology. Volume 57, Issue 3, May 2000, Pages 265-271
    '''
    alleles = snv.qualifiers['alleles']
    read_groups = snv.qualifiers['read_groups']
    alleles_reads = []
    for allele in alleles:
        alleles_reads.append(_allele_count(allele, alleles,
                                           read_groups=read_groups,
                                           group_kind=group_kind,
                                           groups=groups ))

    if len(alleles_reads) == 1:
        pic = 0
    else:
        total_reads = sum(alleles_reads)

        # we need at least 4 reads to calculate pic
        if total_reads < 4:
            pic = None
        else:
            first_element = 0
            second_element = 0
            for index, num_allele in enumerate(alleles_reads):
                first_element_part = ((num_allele*(num_allele - 1))/
                                      (total_reads*(total_reads - 1)))
                first_element += first_element_part

                for num_allele in alleles_reads[index+1:]:
                    second_element_part = (first_element_part*
                                           ((num_allele*(num_allele - 1))/
                                            ((total_reads - 2)*
                                             (total_reads - 3))))
                    second_element += second_element_part
            pic = 1 - first_element - second_element
    return pic

def annotate_pic(snv):
    'It annotates the pic'
    pic = calculate_pic(snv)
    snv.qualifiers['pic'] = pic


def calculate_heterozygosity(snv, ploidy, group_kind=None,
                             groups=None):
    '''It calculates the estimator of heterozygosity, given a list with the
    number of times that each allele has been read.

    heterozygosity(estimator) = (n/(n-1))(1-summatory(xi**2))

    xi = gene frequency of the i-th allele
    n = total number of reads

    Formula taken from "SAMPLING VARIANCES OF HETEROZYGOSITY AND GENETIC
    DISTANCE" by MASATOSHI NEI and A. K. ROYCHOUDHURY.
    Genetics 76: 379-390 February, 1974.

    If the number of individuals is less than 50, formula has to be corrected:

    heterozygosity(estimator) = (2n/(2n-1))(1-summatory(xi**2))

    Taken from "ESTIMATION OF AVERAGE HETEROZYGOSITY AND GENETIC DISTANCE FROM
    A SMALL NUMBER OF INDIVIDUALS" by MASATOSHI NEI.
    Genetics 89 : 583-590 July, 1978.
    '''
    alleles = snv.qualifiers['alleles']
    read_groups = snv.qualifiers['read_groups']
    alleles_reads = []
    for allele in alleles:
        alleles_reads.append(_allele_count(allele, alleles,
                                           read_groups=read_groups,
                                           group_kind=group_kind,
                                           groups=groups))

    if len(alleles_reads) == 1:
        heterozygosity = 0
    else:
        total_reads = sum(alleles_reads)
        if total_reads == 0:
            heterozygosity = None
        else:
            sum_  = 0
            for num_allele in alleles_reads:
                allele_freq = num_allele/total_reads
                sum_ += allele_freq**2

            if total_reads/ploidy < 50:
                heterozygosity = ((2*total_reads)/((2*total_reads) - 1))*(1 - sum_)
            else:
                heterozygosity = ((total_reads)/((total_reads) - 1))*(1 - sum_)
    return heterozygosity

def annotate_heterozygosity(snv, ploidy):
    'It annotates the heterozigosity'
    heterozygosity = calculate_heterozygosity(snv, ploidy)
    snv.qualifiers['heterozygosity'] = heterozygosity
