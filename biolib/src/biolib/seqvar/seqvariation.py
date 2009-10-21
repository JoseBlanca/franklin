'''
Created on 2009 mar 25

@author: peio
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

from biolib.biolib_seqio_utils import temp_fasta_file
from biolib.biolib_cmd_utils import call
from biolib.seqs import SeqWithQuality
from biolib.biolib_seqio_utils import FileSequenceIndex
from biolib.collections_ import item_context_iter

SNP = 0
INSERTION = 1
DELETION = 2
INVARIANT = 3
INDEL = 4
COMPLEX = 5

COMMON_ENZYMES = ['ecori', 'smai', 'bamhi', 'alui', 'bglii',
                  'sali', 'bgli', 'clai', 'bsteii', 'taqi',
                  'psti', 'pvuii', 'hindiii', 'ecorv', 'xbai',
                  'haeiii', 'xhoi', 'kpni', 'scai', 'banii',
                  'hinfi', 'drai', 'apai', 'asp718']

def _allele_reads_compare(allele1, allele2):
    'It returns which allele has been read more times'
    return allele2['reads'] - allele1['reads']

class Snv(object):
    '''
    This class is used to represent a sequence variation respect a reference.
    The variation can be: Snp, insertion, deletion or non_variant.
    '''
    def __init__(self, reference, location, per_lib_info=None, name=None):
        '''It initializes with an alleles dict and a reference.

        The lib_allele is a list of dictionaries. Each of the dictionaries
        contains:
            - library : Library name. Not needed but recommended
            - alleles : alleles should be given with a list.
                The allele info should be a dict with the following info:
                - allele: the sequence variation
                - reads: the number of reads that support them as values.
                - kind: the type of variation
                - qualities: a list of qualities
        The reference is the reference sequence, can be an str or a SeqRecord.
        The location must be an int referencing the position in the reference
        sequence.
        '''
        self.reference = reference
        self.name = name
        self.location = int(location)
        self.annotations = {}
        if per_lib_info is None:
            per_lib_info = []
        self.per_lib_info = per_lib_info

        # we have to sort alleles from each library
        for library in self.per_lib_info:
            library['alleles'] = sorted(library['alleles'],
                                        _allele_reads_compare)

        #for every library we have to add the 'annotations'
        for library in per_lib_info:
            if 'annotations' not in library:
                library['annotations'] = {}

    def _get_kind(self):
        '''It returns the kind of variation.

        The kind of variation depends on the kind of the kinds f each allele
         library.
        '''
        kind = INVARIANT
        for alleles_info in self.per_lib_info:
            alleles = alleles_info['alleles']
            for allele_info in alleles:
                al_kind = allele_info['kind']
                kind = calculate_kind(al_kind, kind)
        return kind
    kind = property(_get_kind)

    def copy(self, per_lib_info=None, reference=None, name=None, location=None,
                                                            annotations=None):
        '''Given a seqvariation in returns a new seqvariation with the new data
         changes'''
        if per_lib_info is None:
            per_lib_info = self.per_lib_info
        if reference is None:
            reference = self.reference
        if name is None:
            name = self.name
        if location is None:
            location = self.location

        snv = self.__class__(per_lib_info=per_lib_info, reference=reference,
                             name=name, location=location)
        if annotations is None:
            snv.annotations = self.annotations
        else:
            snv.annotations = annotations
        return snv

    def __str__(self):
        'It print some minimal info'
        to_print = '%s: %d' % (self.reference, self.location)
        return to_print

    def __repr__(self):
        'It prints an evaluable representation'
        to_print  = '%s(\nreference=%s, location=%s,\n' % \
            (self.__class__.__name__, repr(self.reference),repr(self.location))

        to_print += '\tper_lib_info=[\n'
        for alleles_in_a_lib in self.per_lib_info:
            to_print += '\t\t{'
            for key, value in alleles_in_a_lib.items():
                if key != 'alleles':
                    to_print += '\t\t"%s": %s,' % (key, repr(value))
                else:
                    to_print += '"alleles" :[\n'
                    for allele in value:
                        to_print += '\t\t\t%s,\n' % repr(allele)
                    to_print += '\t\t\t],\n'
            to_print += '\t\t},\n'
        to_print += '])\n\n'
        return to_print


def cap_enzymes(snv, all_enzymes=False):
    '''Given an svn it returns the list of restriction enzymes that distinguish
    between their alleles.'''
    if 'cap_enzymes' in snv.annotations:
        return snv.annotations

    #which alleles do we have?
    alleles = set()
    for lib in snv.per_lib_info:
        for allele in lib['alleles']:
            alleles.add(repr((allele['allele'], allele['kind'])))
    #for every pair of different alleles we have to look for differences in
    #their restriction maps
    enzymes = set()
    alleles = list(alleles)
    reference = snv.reference
    location = snv.location
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

    enzymes = list(enzymes)
    snv.annotations['cap_enzymes'] = enzymes
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
    ref = reference
    loc = location
    def create_sequence(name, allele, kind):
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
        seq = SeqWithQuality(name=name, seq=seq)
        return seq
    seq1 = create_sequence('seq1', allele1, kind1)
    seq2 = create_sequence('seq2', allele2, kind2)

    enzymes1 = _remap_run(seq1, all_enzymes)
    enzymes2 = _remap_run(seq2, all_enzymes)

    enzymes = set(enzymes1).symmetric_difference(set(enzymes2))

    return enzymes

def _remap_run(seq, all_enzymes):
    '''this command runs remap EMBOSS binary and returns ...'''
    # Minimun length of the restriction enzyme recognition site
    sitelen = 4
    seq_file = temp_fasta_file(seq)
    seq_filename = seq_file.name

    if all_enzymes:
        enzymes = 'all'
    else:
        enzymes = ",".join(COMMON_ENZYMES)

    cmd = ['remap', '-sequence', seq_filename, '-enzymes', enzymes,
           '-sitelen' , str(sitelen), 'stdout']

    try:
        stdout, stderr, retcode  = call(cmd)
    except OSError:
        raise OSError('remap binary does not exits or it is not in path')

    if retcode:
        raise RuntimeError('remap error: '+ stderr)

    return  _parse_remap_output(stdout)

def _parse_remap_output(remap_output):
    ''' It takes the remap output and it returns a set list with the enzymes
     that cut there'''
    section = ''
    enzymes = []
    for line in remap_output.split('\n'):
        line = line.strip()
        if line.isspace() or len(line) < 2:
            continue
        if section == 'cut':
            if line.startswith('#'):
                section = ''
            else:
                enzymes.append(line.split()[0])

        if line.startswith('# Enzymes that cut'):
            section = 'cut'
            continue

    return enzymes

def reference_variability(snv, context, window=None):
    'It calculates the variability of thh reference of the snv'
    #how many snps are in the window?
    snv_quantity = len(context)
    if window is None:
        if 'seq' not in dir(snv.reference):
            raise ValueError('The reference should be a seqRecord')
        window = len(snv.reference)
    return snv_quantity / float(window) * 100.0

def _maf_for_alleles_in_lib(alleles):
    'It returns the maf for the given alleles'
    most_freq_reads = alleles[0]['reads']
    tot_reads = 0
    for allele in alleles:
        tot_reads += allele['reads']
    maf = most_freq_reads / float(tot_reads)
    return maf

#functions to characterize the sequence variations
def _pic_for_alleles_in_lib(alleles):
    '''It calculates and returns the Polymorphic Information Content.

    The PIC was defined in Botstein 1980. Am. J. Hum. Genet. 32, 314 331 as  the
    probability that a given marker genotype of an offspring of an affected
    parent will allow deduction of the parental genotype at the marker locus.
    The calculation is done following Shete et al. 2000 Theoretical Population
    Biology 57, 265 271.
    '''
    #pylint: disable-msg=C0103
    def _pic_sum_1(alleles, num_reads, num_alleles):
        '''It returns the first summation for the pic calculation'''
        # P is the frecuency one allele have been read
        suma  = 0
        for i in range(num_alleles):
            frec = (alleles[i] / num_reads) ** 2
            suma += frec
        return suma

    def _pic_sum_2(alleles, num_reads, num_alleles):
        '''It returns the second summation for the pic calculation '''
        suma = 0
        for i in range(num_alleles - 1):
            freci = (alleles[i] / num_reads) ** 2
            for j in range(i + 1, num_alleles):
                frecj = (alleles[j] / num_reads) ** 2
                suma += freci * frecj
        return suma

    #the alleles can have the count or a list with the alleles, we make sure
    #that all have a count, and we convert the dict to a list
    alleles = [allele['reads'] for allele in alleles]
    #how many reads are in total?
    num_reads = float(sum(alleles))
    #how many alleles are in total
    num_alleles = len(alleles)
    sum_1 = _pic_sum_1(alleles, num_reads, num_alleles)
    sum_2 = _pic_sum_2(alleles, num_reads, num_alleles)

    return 1.0 - sum_1 - ( 2 * sum_2)

SVN_ANNOTATION_CALCULATORS = {'maf': _maf_for_alleles_in_lib,
                              'pic': _pic_for_alleles_in_lib}

def _calculate_libs_annotation(snv, kind):
    'It calculates an annotation for every library in the snv'
    annots = []
    for lib in snv.per_lib_info:
        if kind in lib['annotations']:
            value = lib['annotations'][kind]
        else:
            alleles = lib['alleles']
            value = SVN_ANNOTATION_CALCULATORS[kind](alleles)
        annots.append(value)
        #we store the value in the annotations
        lib['annotations'][kind] = value
    return annots

def major_allele_frequency(snv):
    'It returns a list with the frequencies of the maf in all libraries'
    return _calculate_libs_annotation(snv, 'maf')

def pic(snv):
    'It returns a list with the pics in all libraries'
    return _calculate_libs_annotation(snv, 'pic')


def _snvs_in_file(snv_fhand):
    'It reads the snv evalable file and it yields snvs'
    snv_buffer = ''
    for line in snv_fhand:
        line = line.strip()
        if not line and snv_buffer:
            snv = eval(snv_buffer)
            yield snv
            snv_buffer = ''
        snv_buffer += line
def snvs_in_file(snv_fhand, ref_fhand=None):
    'It reads the snv evalable file and it yields snvs'
    snvs = _snvs_in_file(snv_fhand)
    if ref_fhand:
        snvs = add_reference_to_svns(snvs, ref_fhand)
    return snvs

def add_reference_to_svns(snvs, ref_fhand):
    'It adds the reference to the svn, it yields the complete svn'
    references_index = FileSequenceIndex(ref_fhand)
    for snv in snvs:
        snv.reference = references_index[snv.reference]
        yield snv

def svn_contexts_in_file(snv_fhand, ref_fhand=None):
    'It reads an svn file and it yields (svn, context) tuples'
    return item_context_iter(snvs_in_file(snv_fhand, ref_fhand))

def calculate_kind(kind1, kind2):
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

