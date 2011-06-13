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

from collections import defaultdict
from copy import copy
import math

try:
    import pysam
except ImportError:
    pass

from Bio.SeqFeature import FeatureLocation
from Bio.Restriction import Analysis, CommOnly, RestrictionBatch
from franklin.seq.seqs import SeqFeature, get_seq_name
from franklin.utils.misc_utils import get_fhand
from franklin.sam import create_bam_index, get_read_group_info

DELETION_ALLELE = '-'
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

def _get_raw_allele_from_read(aligned_read, index):
    'It returns allele, quality, is_reverse'
    allele = aligned_read.seq[index].upper()
    if aligned_read.qual:
        qual = _qualities_to_phred(aligned_read.qual[index])
    else:
        qual = None
    return allele, qual

def _get_segments_from_cigar(ref_pos, cigar, read_len):
    '''It returns two lists (reference and read) in which the firsts nucleotides
     of the different cigar categories are given.

     67890  12345
     ATCGA--GATCG
     atcGATCG--CG
        12345  67

    CIGAR = 3H2M2I1M2D2M

    ref_segments = [9, None, 11, 12, 14]
    read_segments = [1, 3, 5, None, 6]

    It also returns the limits of the aligned reference.

    It also returns a list with the cigar category for each segment.
    '''

    #We ignore hard clipped nucleotides ('H')
    cigar_elements = []
    for element in range(len(cigar)):
        if cigar[element][0] != HARD_CLIP:
            cigar_elements.append(cigar[element])

    ref_segments = []
    read_segments = []
    ref_start = ref_pos
    segment_type = []
    segment_lens = []

    read_pos = 0

    for element in cigar_elements:
        if element[0] == SOFT_CLIP:
            read_pos += element[1]
            continue
        elif element[0] == MATCH:
            segment_type.append(MATCH)
            ref_segments.append(ref_pos)
            read_segments.append(read_pos)
            segment_lens.append(element[1])
            ref_pos += element[1]
            read_pos += element[1]
            ref_start += element[1]
            ref_end = ref_start - 1
        elif element[0] == INSERTION:
            segment_type.append(INSERTION)
            ref_segments.append(None)
            read_segments.append(read_pos)
            segment_lens.append(element[1])
            read_pos += element[1]
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
            ref_start += element[1]
            ref_end = ref_start - 1

    ref_start = ref_segments[0]
    ref_limits = [ref_start, ref_end]
    return (ref_segments, read_segments, sorted(ref_limits),
            segment_type, segment_lens)

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

def _get_insertion(segment_index, segment_type, read_pos, pileup_read,
                  aligned_read):

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

        indel_length = pileup_read.indel
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

def _get_alleles_from_read(ref_allele, ref_pos, pileup_read):
    '''It returns an allele from the read.

    It returns a list with the alleles in the given position.
    The returned allele can be an empty list if we're in a deletion.
    If the position holds an insertion it will return two alleles, the
    insertion and the nucleotide at that position.
    '''

    alleles = []
    aligned_read = pileup_read.alignment
    begin_pos_read_in_ref = aligned_read.pos
    read_len = len(aligned_read.seq)
    cigar = aligned_read.cigar

    (ref_segments, read_segments, ref_limits, segment_types,
    segment_lens) = _get_segments_from_cigar(begin_pos_read_in_ref, cigar,
                                             read_len)
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
                                               next_read_pos, pileup_read,
                                               aligned_read)
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
                qual0 = _qualities_to_phred(qual0)
                qual1 = aligned_read.qual[read_pos2]
                qual1 = _qualities_to_phred(qual1)
                qual = min((qual0, qual1))
            else:
                qual = None
            kind = DELETION
            alleles.append((allele, kind, qual, is_reverse))
        if segment_pos ==IN_FIRST_AND_LAST or segment_pos == IN_LAST_POS:
            #Is there an insertion in the next position?
            allele, kind, qual = _get_insertion(segment_index, segment_types,
                                               read_pos1, pileup_read,
                                               aligned_read)
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
                                           read_pos1, pileup_read,
                                           aligned_read)
        if kind is not None:
            alleles.append((allele, kind, qual, is_reverse))

    elif segment_types[segment_index] == INSERTION:
        pass    #if we're in an insertion, it is returned in the last position
                #of the previous match segment
    return alleles
def _qualities_to_phred(quality):
    'It transforms a qual chrs into a phred quality'
    if quality is None:
        return None
    phred_qual = []
    for char in quality:
        phred_qual.append(ord(char) - 33)
    if quality[0] == 93:  #the character used for unknown qualities
        phred_qual = None
    else:
        phred_qual = sum(phred_qual) / len(phred_qual)
    return phred_qual

def _add_allele(alleles, allele, kind, read_name, read_group, is_reverse, qual,
                mapping_quality, readgroup_info):
    'It adds one allele to the alleles dict'
    key = (allele, kind)
    if key not in alleles:
        alleles[key] = {'read_groups':[], 'orientations':[],
                        'qualities':[], 'mapping_qualities':[]}
    allele_info = alleles[key]
    allele_info['read_groups'].append(read_group)
    allele_info['orientations'].append(not(is_reverse))
    allele_info['qualities'].append(qual)
    allele_info['mapping_qualities'].append(mapping_quality)

def _normalize_read_edge_conf(read_edge_conf):
    'It returns a dict with all valid keys'
    platforms = ('454', 'sanger', 'illumina')
    if read_edge_conf is None:
        read_edge_conf = {}
    for platform in platforms:
        if platform not in read_edge_conf:
            read_edge_conf[platform] = (None, None)
    return read_edge_conf

def _snvs_in_bam(bam, reference, min_quality, default_sanger_quality,
                 min_mapq, min_num_alleles, max_maf, min_num_reads_for_allele,
                 read_edge_conf=None, default_bam_platform=None):
    'It yields the snv information for every snv in the given reference'

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
            if read_edge_conf and platform in read_edge_conf:
                edge_left, edge_right = read_edge_conf[platform]

                #if we're in the edge region to be ignored we continue to the next
                #read, because there's no allele to add for this one.
                if ((edge_left  is not None and edge_left >= read_pos) or
                    (edge_right is not None and edge_right <= read_pos)):
                    continue

            alleles_here = _get_alleles_from_read(ref_allele, ref_pos,
                                                  pileup_read)
            for allele in alleles_here:
                allele, kind, qual, is_reverse = allele
                _add_allele(alleles, allele, kind, read_name, read_group,
                    is_reverse, qual, read_mapping_qual,
                    read_groups_info)

        #remove N
        _remove_alleles_n(alleles)

        #add default sanger qualities to the sanger reads with no quality
        _add_default_sanger_quality(alleles, default_sanger_quality,
                                    read_groups_info)

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
    quals.sort(lambda x, y: int(y - x))

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

    #we remove from the read_groups the ones not used in this snv
    new_read_groups = {}
    for read_group, info in snv['read_groups'].items():
        if read_group in used_read_groups:
            new_read_groups[read_group] = info

    snv['read_groups'] = new_read_groups

    return snv

def create_snv_annotator(bam_fhand, min_quality=45, default_sanger_quality=25,
                         min_mapq=15, min_num_alleles=1, max_maf=0.9,
                         read_edge_conf=None, default_bam_platform=None,
                         min_num_reads_for_allele=2):
    'It creates an annotator capable of annotating the snvs in a SeqRecord'

    #the bam should have an index, does the index exists?
    bam_fhand = get_fhand(bam_fhand)
    create_bam_index(bam_fpath=bam_fhand.name)
    read_edge_conf = _normalize_read_edge_conf(read_edge_conf)

    bam = pysam.Samfile(bam_fhand.name, 'rb')

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
            feat = SeqFeature(location=FeatureLocation(location, location),
                              type=type_,
                              qualifiers=qualifiers)
            sequence.features.append(feat)
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

def snvs_in_window(snv, snvs, window, snv_type=None):
    'it gets all the snvs  in a window taking a snv as reference'
    num_of_snvs = 0
    location = int(str(snv.location.start))
    left_margin = location - (window / 2)
    rigth_margin = location + (window / 2)
    for snv in snvs:

        location = int(str(snv.location.start))
        if location > left_margin and location < rigth_margin:
            if  snv_type is None:
                num_of_snvs += 1
            else:
                type_ = calculate_snv_kind(snv)
                if ((snv_type == type_) or
                    (snv_type == INDEL and type_ in(INSERTION, DELETION))):
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
                           min_reads_per_allele=None):
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
                                       maf=maf)
        variable_per_read_group.append(variable)

    if in_all_groups:
        return all(variable_per_read_group)
    else:
        return any(variable_per_read_group)

def invariant_in_groupping(snv, group_kind, groups, in_union=True,
                           in_all_groups=True, reference_free=True):
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
                                         reference_free=reference_free)
        invariable_per_read_group.append(invariable)

    if in_all_groups:
        return all(invariable_per_read_group)
    else:
        return any(invariable_per_read_group)


#def _is_rg_invariable(alleles, reference_free=True):
#    'It checks if the allele is variable'
#    allele_count = len(alleles)
#
#    ref_in_alleles = False
#    for allele in alleles.keys():
#        if allele[1] == INVARIANT:
#            ref_in_alleles = True
#
#    if allele_count == 1:
#        if reference_free:
#            return True
#        elif not reference_free and ref_in_alleles:
#            return True
#
#    return False

def _is_rg_variable(alleles, reference_free=True, maf=None):
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

    if maf and maf < calc_maf(alleles):
        return False

    return True

def calc_maf(alleles):
    'It calculates the maf'
    values = alleles.values()
    major = max(values)
    sum_  = sum(values)
    return  major/float(sum_)

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



def variable_in_groupping_old(group_kind, feature, groups, in_union=False,
                           in_all_groups=True):
    'It looks if the given snv is variable for the given groups'

    alleles = _get_alleles_for_group(feature.qualifiers['alleles'],
                                     groups, group_kind,
                                     feature.qualifiers['read_groups'])
    if not alleles:
        return None

    if in_union:
        alleles = _aggregate_alleles(alleles)

    variable_in_read_groups_ = []
    for allele_list in alleles.values():
        variable_in_read_groups_.append(True if len(allele_list) > 1 else False)

#    #For the case in which there are no alleles
#    if not variable_in_read_groups_:
#        return None

    if in_all_groups:
        return all(variable_in_read_groups_)
    else:
        return any(variable_in_read_groups_)

def _aggregate_alleles_old(alleles):
    'It joins all alleles for the read groups into one'
    aggregate = set()
    for allele_list in alleles.values():
        aggregate = aggregate.union(allele_list)
    return {None: aggregate}

def _get_alleles_for_group_old(alleles, groups, group_kind='read_groups',
                           read_groups=None):
    '''It gets the alleles from the given items of type:key, separated by items.
    For example, if you give key rg and items rg1, rg2, it will return
    alleles separated in rg1 and rg2 '''

    alleles_for_groups = {}
    for allele, alleles_info in alleles.items():
        for read_group in alleles_info['read_groups']:
            group = _get_group(read_group, group_kind, read_groups)
            if group not in groups:
                continue
            if not group in alleles_for_groups:
                alleles_for_groups[group] = set()
            alleles_for_groups[group].add(allele)
    return alleles_for_groups
