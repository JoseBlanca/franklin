
'''
Created on 22/09/2009

@author: jose
'''

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

from biolib.biolib_seqio_utils import FileSequenceIndex
from biolib.collections_ import item_context_iter
from biolib.seqvar.seqvariation import (SNP, INSERTION, DELETION,
                                        INVARIANT, Snv)
import copy

def _get_alleles_from_line(line_split):
    'It gets allele from each line of the pileup'
    cromosome, position, ref_base, coverage, read_bases, qual = line_split
    # get alleles from the string
    try:
        alleles, qualities = _get_alleles(ref_base, read_bases, qual)
    except RuntimeError:
        msg = 'malformed line in pileup: %s %s\n' % (cromosome, position)
        raise ValueError(msg)

    alleles, qualities, orientations = _group_alleles(alleles, qualities)
    alleles = _get_allele_type(ref_base, alleles, qualities, orientations)
    return alleles

def _get_allele_type(ref_base, alleles, qualities, orientations):
    '''It gets the type (snp, deletion, insertion) for each allele and removes
    the + and -'''
    new_alleles = []
    for allele, num_reads in alleles.items():
        allele_info = {}
        allele_info['reads']   = num_reads
        allele_info['qualities'] =  qualities[allele]
        allele_info['orientations'] = orientations[allele]
        if allele.startswith('+'):
            allele_info['allele'] = allele[1:]
            allele_info['kind']   = INSERTION
        elif allele.startswith('-'):
            allele_info['allele'] = allele[1:]
            allele_info['kind']   = DELETION
        elif allele == ref_base:
            allele_info['allele'] = allele
            allele_info['kind']   = INVARIANT
        else:
            allele_info['allele'] = allele
            allele_info['kind']   = SNP
        new_alleles.append(allele_info)
    return new_alleles


def _is_seq_var(ref_base, alleles, min_reads):
    '''This looks if the given data is a seq variations, it return false if it
     is not a seq var and return allele, qual_allele if it is a seq_var'''

    min_alleles = 2 # At least it needs a
    allele_nucleotides = []
    for allele in alleles:
        allele_nucleotides.append(allele['allele'])
    if len(alleles) < min_alleles and ref_base in allele_nucleotides:
        return False
    return True

def _group_alleles(alleles, quals):
    '''It converts the list of seqs and quals into dicts with the allele as key.
    allele_group group alleles and counts them and quality_group gives
    the quality for each allele'''
    alleles_dict = {}
    qual_dict = {}
    orientations_dict = {}
    for index, allele in enumerate(alleles):
        al_orig = allele
        allele = allele.upper()
        if allele not in alleles_dict:
            alleles_dict[allele] = 0
            qual_dict[allele] = []
            orientations_dict[allele] = []
        if al_orig == allele:
            orientations_dict[allele].append(True)
        else:
            orientations_dict[allele].append(False)
        alleles_dict[allele] += 1
        qual_dict[allele].append(quals[index])
    return alleles_dict, qual_dict, orientations_dict

def _qualities_to_phred(qualities):
    'It transforms a list of chrs into a phred list'
    new_qualities = []
    for qual in qualities:
        if qual is not None:
            qual = ord(qual) - 33
        if qual == 93:  #the character used for unknown qualities
            qual = None
        new_qualities.append(qual)
    return new_qualities

def _get_alleles(ref_base, alleles, qualities):
    '''Given the sequence and qualities strings it returns two lists with the
    sequence and quality'''
    alleles = list(alleles)
    qualities = list(qualities)
    def indel_span(pos):
        'how many positions the indel definitions takes'
        number = alleles[pos]
        while True:
            pos += 1
            try:
                new_number = alleles[pos]
            except IndexError:
                raise RuntimeError('Malformed pileup line')
            if not new_number.isdigit():
                break
            number += new_number
        return int(number)
    al_pos = 0
    qual_pos = 0
    pos_delta = 1
    qual_delta = 1
    filtered_alleles = []
    filtered_qualities = []
    ignore_qual = False
    while al_pos < len(alleles):
        item = alleles[al_pos:al_pos + pos_delta]
        if item == ['.']:
            item = [ref_base]
        elif item == [',']:
            item = [ref_base.lower()]
        if item == ['$']:
            item = []
            qual_delta = 0
            ignore_qual = True
        elif item == ['*']:
            item = []
            ignore_qual = True  #the deletion has a quality associated we ignore
        elif item == ['^']:
            pos_delta = 2
            qual_delta = 0
            item = []
            ignore_qual = True
        elif item == ['+'] or item == ['-']:
            pos_delta = indel_span(al_pos + 1)
            al_pos += len(str(pos_delta)) + 1
            item.extend(alleles[al_pos:al_pos + pos_delta])
            qual_delta = 0 #no quality for the indels
        if item:
            filtered_alleles.append(''.join(item))
        al_pos += pos_delta
        pos_delta = 1
        if not ignore_qual and qual_delta:
            qual = qualities[qual_pos]
            filtered_qualities.append(qual)
        elif qual_delta == 0 and not ignore_qual:
            filtered_qualities.append(None)
        qual_pos += qual_delta
        qual_delta = 1
        ignore_qual = False
    filtered_qualities = _qualities_to_phred(filtered_qualities)
    return filtered_alleles, filtered_qualities

def save_seqvars_positions(seq_vars, outfhand):
    '''It takes a seqvar iterator and writes the required position file.
     each line of this file is a tuple separated by a tab of reference and
     postion'''
    for seq_var in seq_vars:
        location = seq_var[0].location
        try:
            name = seq_var[0].reference.name
        except AttributeError:
            name =  seq_var[0].reference
        outfhand.write('%s\t%d\n' % (name, location))
    outfhand.close()


def _fill_cache(pileups, cache):
    'It replenishes the cache with the cache'
    for index in range(len(pileups)):
        if not cache[index]:
            #if this pileup cache is fill, nothing to do
            #else we get one line
            try:
                line = pileups[index].next().strip()
            except StopIteration:
                line = ''
            cache[index] = line.split()

def _index_lower(index1, index2):
    'Compares two Snv ref, locations and returns True if the 1 is lower than 2'
    if (index1[0] < index2[0] or
        (index1[0] == index2[0] and index1[1] < index2[1])):
        return True
    return False

def _lowest_location(pileup_lines):
    'Given a list with pileup lines it returns the lowest one'
    lowest_index = (None, None)
    for pileup_line in pileup_lines:
        this_index = pileup_line[:2]
        if not this_index:  #an empty index when the pileup file is over
            continue
        if (lowest_index == (None, None) or
            not _index_lower(lowest_index, this_index)):
            lowest_index = this_index
    return lowest_index

def _pop_lines_from_cache(index, cache):
    'It returns the lines from the cache that have the given index'
    lines = []
    all_empty = True
    for cache_index in range(len(cache)):
        if cache[cache_index][:2] == index:
            lines.append(cache[cache_index])
            cache[cache_index] = []
            all_empty = False
        else:
            lines.append([])
    if all_empty:
        raise StopIteration
    return lines

def _locations_in_pileups(pileups):
    'It yields the equivalent positions from different pileup files'
    cache = [[]] * len(pileups)
    _fill_cache(pileups, cache)
    while True:
        lowest_index = _lowest_location(cache)
        lines_in_index = _pop_lines_from_cache(lowest_index, cache)
        yield lines_in_index
        _fill_cache(pileups, cache)


def _alleles_in_pileups(line_in_pileups):
    'It gets the alleles of each of the pileups.'
    lib_alleles = []
    for line_split in line_in_pileups:
        if line_split:
            try:
                alleles = _get_alleles_from_line(line_split)
            except ValueError:
                #TODO Add a log here
                alleles = []
        else:
            alleles = []
        lib_alleles.append(alleles)
    return lib_alleles

def agregate_alleles(alleles):
    'It aggreates the alleles from a list of alleles (one for every lib)'
    ag_alleles = {}
    for lib_alleles in alleles:
        for allele in lib_alleles:
            base = allele['allele']
            kind = allele['kind']
            if (base, kind) not in ag_alleles:
                ag_alleles[(base, kind)] = copy.deepcopy(allele)
            else:
                ag_al = ag_alleles[(base, kind)]
                ag_al['reads'] += allele['reads']
                ag_al['qualities'].extend(allele['qualities'])
    return ag_alleles.values()

def snvs_in_sam_pileups(pileups, libraries, references=None, min_num=None):
    '''It yields Snvs from the given pileups'''
    if references is not None:
        references_index = FileSequenceIndex(references)

    # check if the pileups are well formed
    for pileup in pileups:
        if not _check_pileup(pileup):
            raise ValueError(' malformed pile up file: %s' % pileup.name)

    for lines_in_pileups in _locations_in_pileups(pileups):
        lib_alleles = _alleles_in_pileups(lines_in_pileups)
        ag_alleles = agregate_alleles(lib_alleles)
        #which is the ref_base
        non_empty_line = None
        for line_in_pileup in lines_in_pileups:
            if line_in_pileup:
                non_empty_line = line_in_pileup
                break
        ref_base = non_empty_line[2]
        if not _is_seq_var(ref_base, ag_alleles, min_num):
            continue
        #now we can build the Snv
        if references is not None:
            reference = references_index[non_empty_line[0]]
        else:
            reference = non_empty_line[0]
        location = non_empty_line[1]
        per_lib_info = []
        for lib_index, lib_allele in enumerate(lib_alleles):
            library = libraries[lib_index]
            if lib_allele:
                per_lib_info.append({'alleles':lib_allele,
                                     'library':library,
                                     'annotations':{},
                                     },
                                    )
        yield Snv(reference=reference, location=location,
                  per_lib_info=per_lib_info)

def snv_contexts_in_sam_pileup(pileups, libraries, min_num=None, window=None,
                       references=None):
    '''This function takes the seqvar iterator of the sam pileup, and it return
    an iterator with a tuple of seqvar and its context'''
    seqvar_iter = snvs_in_sam_pileups(pileups, libraries=libraries,
                                       min_num=min_num, references=references)
    seq_var_with_context_iter = item_context_iter(seqvar_iter, window=window)
    for seq_var_contex in seq_var_with_context_iter:
        yield seq_var_contex

def _check_pileup(pileup):
    'It check if the pileup is well formed and its data is correct'
    notdotcomma = 0
    total_nt    = 0
    total_lines = 0
    ref_base_n  = 0
    for i, line in enumerate(pileup):
        if i >= 100:
            break
        items = line.split()
        if items[2].upper() == 'N':
            ref_base_n += 1
        total_lines += 1
        for nucl in items[4]:
            if nucl not in [',', '.', '$', '^']:
                notdotcomma += 1
            total_nt += 1
    pileup.seek(0)
    if ((notdotcomma / float(total_nt) >= 0.5) or
        (ref_base_n / float(total_lines) >= 0.3)):
        return False
    return True

