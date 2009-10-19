
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

import sys

from biolib.biolib_seqio_utils import FileSequenceIndex
from biolib.collections_ import item_context_iter
from biolib.seqvar.seqvariation import (SNP, INSERTION, DELETION,
                                        INVARIANT, Snv)

def _seqvars_in_sam_pileup(pileup, min_num=None, required_positions=None,
                           references=None, library=None):
    '''This function takes from a sam pileup format file all the position that
    are interesting'''
    if references is not None:
        references_index = FileSequenceIndex(references)
    if min_num == None:
        min_num = 2
    for line in pileup:
        if line.isspace():
            continue

        cromosome, position, ref_base, coverage, read_bases, qual = line.split()

        if required_positions and not required_positions[cromosome, position]:
            continue

        # get alleles from the string
        try:
            alleles, qualities = _get_alleles(ref_base, read_bases, qual)
        except RuntimeError:
            sys.stderr.write('malformed line in pileup: %s %s\n' % (cromosome,
                                                                    position))
            continue
        alleles, qual_grouped = _group_alleles(alleles, qualities)
        alleles = _get_allele_type(ref_base, alleles, qual_grouped)

        if required_positions:
            yield_snv = required_positions[cromosome, position]
        else:
            yield_snv = _is_seq_var(coverage, ref_base, alleles, min_num)
        if yield_snv:
            if references is not None:
                cromosome = references_index[cromosome]
            lib_alleles = {}
            lib_alleles['alleles'] = alleles
            if library is not None:
                lib_alleles['library'] = library
            #print lib_alleles
            yield Snv(per_lib_info=[lib_alleles],
                      location=int(position) - 1, reference=cromosome)

def seqvars_in_sam_pileup(pileup, min_num=None, window=None, library=None,
                          required_positions=None, references=None):
    '''This function takes the seqvar iterator of the sam pileup, and it return
    an iterator with a tuple of seqvar and its context'''
    seqvar_iter = _seqvars_in_sam_pileup(pileup, min_num, required_positions,
                                         references, library)
    seq_var_with_context_iter = item_context_iter(seqvar_iter, window=window)
    for seq_var_contex in seq_var_with_context_iter:
        yield seq_var_contex


def _get_allele_type(ref_base, alleles, qual_grouped):
    '''It gets the type (snp, deletion, insertion) for each allele and removes
    the + and -'''
    new_alleles = []
    for allele, num_reads in alleles.items():
        allele_info = {}
        allele_info['reads']   = num_reads
        allele_info['qualities'] =  qual_grouped[allele]
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


def _is_seq_var(coverage, ref_base, alleles, min_reads):
    '''This looks if the given data is a seq variations, it return false if it
     is not a seq var and return allele, qual_allele if it is a seq_var'''

    min_alleles = 2 # At least it needs a
    if coverage < min_reads:
        return False
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
    for index, allele in enumerate(alleles):
        allele = allele.upper()
        if allele not in alleles_dict:
            alleles_dict[allele] = 0
            qual_dict[allele] = []
        alleles_dict[allele] += 1
        qual_dict[allele].append(quals[index])
    return alleles_dict, qual_dict

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
        if item == ['.'] or item == [',']:
            item = [ref_base]
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







