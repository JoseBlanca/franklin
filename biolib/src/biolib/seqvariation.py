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
from biolib.biolib_utils import get_start_end
from biolib.biolib_cmd_utils import call
from biolib.seqs import SeqWithQuality

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

class SeqVariation(object):
    '''
    This class is used to represent a sequence variation respect a reference.
    The variation can be: Snp, insertion, deletion or non_variant.
    '''
    def __init__(self, alleles, reference, name=None, location=None):
        '''It initializes with an alleles dict and a reference.

        The alleles should be given with a list.
        The allele info should be a dict with the following info:
            - allele: the sequence variation
            - reads: the number of reads that support them as values.
            - kind: the type of variation
            - qualities: a list of qualities
        The reference is the reference sequence, can be an str or a SeqRecord.
        The location should be an int refering to the location in the reference
        sequence.
        '''
        self.alleles = sorted(alleles, _allele_reads_compare)
        self.reference = reference
        self.name = name
        self.location = location
        self.annotations = {}

    def _get_kind(self):
        '''It returns the kind of variation.

        The kind of variation depends on the kind of alleles.
        '''
        kind = INVARIANT
        for allele_info in self.alleles:
            al_kind = allele_info['kind']
            if kind == INVARIANT and al_kind == DELETION:
                kind = DELETION
            elif kind == INVARIANT and al_kind == INSERTION:
                kind = INSERTION
            elif (kind == DELETION and al_kind == INSERTION or
                 kind == INSERTION and al_kind == DELETION):
                kind = INDEL
            elif ((kind in (DELETION, INSERTION, INDEL) and
                   al_kind == SNP) or
                   (kind == SNP  and al_kind in (INSERTION, INDEL))):
                kind = COMPLEX
            elif kind == INVARIANT and al_kind == SNP:
                kind = SNP
        return kind
    kind = property(_get_kind)

    def __str__(self):
        return '%s %d: %s:' % (self.reference, self.location, str(self.alleles))

    def copy(self, alleles=None, reference=None, name=None, location=None,
                                                            annotations=None):
        '''Given a seqvariation in returns a new seqvariation with the new data
         changes'''
        if alleles is None:
            alleles = self.alleles
        if reference is None:
            reference = self.reference
        if name is None:
            name = self.name
        if location is None:
            location = self.location

        seqvar = self.__class__(alleles=alleles, reference=reference, name=name,
                              location=location)
        if annotations is None:
            seqvar.annotations = self.annotations
        else:
            seqvar.annotations = annotations
        return seqvar

def seqvar_summary(seqvar):
    '''It takes a seqvar and summarizes, returning a library file like format
     string'''
    try:
        contig_name = seqvar.reference.name
    except:
        if len(seqvar.reference) < 30:
            contig_name = seqvar.reference
        else:
            contig_name = ''
    sorted_alleles = allelesseqvar.sorted_alleles()
    second_alleles = allele_count(sorted_alleles[1][1])

    kind               = seqvar.kind()
    loc_start, loc_end = get_start_end(seqvar.location)
    id_snp             = contig_name + '_' +  str(loc_start)
    if 'pic' not in seqvar.annotations:
        calculate_pic(seqvar)
    if 'enzyme' not in seqvar.annotations:
        cap_enzime(seqvar)
    pic                = seqvar.annotations['pic']
    enzymes            = seqvar.annotations['enzymes']

    toprint  = "snp\n"
    toprint += "\tname:%s\n" % id_snp
    toprint += "\tcontig:%s\n" % contig_name
    toprint += "\tstart:%d\n" % loc_start
    toprint += "\tend:%d\n" % loc_end
    toprint += "\tkind:%s\n" % kind
    toprint += "\tannotations: pic:%f\nsecond_allele:%d\nenzymes:%s\n" % \
                                                (pic, second_alleles, enzymes)
    alleles_toprint = ''
    for allele, reads in seqvar.alleles.items():
        reads = [seqvar.alignment[read].name for read in reads]
        alleles_toprint += "%s:%s" % (allele, ','.join(reads))
    toprint += "\talleles:%s\n" % alleles_toprint

    return toprint

#functions to characterize the sequence variations
def calculate_pic(seq_variation):
    '''It calculates and returns the Polymorphic Information Content.

    The PIC was defined in Botstein 1980. Am. J. Hum. Genet. 32, 314 331 as  the
    probability that a given marker genotype of an offspring of an affected
    parent will allow deduction of the parental genotype at the marker locus.
    The calculation is done following Shete et al. 2000 Theoretical Population
    Biology 57, 265 271.
    '''
    # If it is already calculated it returns it
    if 'pic' in seq_variation.annotations:
        return seq_variation.annotations['pic']

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


    alleles = seq_variation.alleles
    #the alleles can have the count or a list with the alleles, we make sure
    #that all have a count, and we convert the dict to a list
    alleles = [allele['reads'] for allele in alleles]
    #how many reads are in total?
    num_reads = float(sum(alleles))
    #how many alleles are in total
    num_alleles = len(alleles)
    sum_1 = _pic_sum_1(alleles, num_reads, num_alleles)
    sum_2 = _pic_sum_2(alleles, num_reads, num_alleles)


    pic = 1.0 - sum_1 - ( 2 * sum_2)
    # We add the pik to the snp annotations
    seq_variation.annotations['pic'] = pic
    return pic

def cap_enzime(snp, all_enzymes=False):
    ''' It looks in the 2 most frecuent alleles if there is each of the enzimes
    cut diferently'''
    # If it is already calculated we return them
    if 'enzymes' in snp.annotations:
        return snp.annotations['enzymes']

    # we look for a CAP between the two most abundant alleles.
    allele1 = snp.alleles[0]
    allele2 = snp.alleles[1]
    kind1 = allele1['kind']
    kind2 = allele2['kind']
    allele1 = allele1['allele']
    allele2 = allele2['allele']

    #we have to build the two sequences
    ref = snp.reference
    loc = snp.location
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

    # We add the enzymes to the snp annotations
    snp.annotations['enzymes'] = enzymes
    return list(enzymes)


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
