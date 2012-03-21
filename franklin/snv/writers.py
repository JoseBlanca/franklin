'''
Created on 19/02/2010

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

from tempfile import NamedTemporaryFile
import datetime, math

from franklin.snv.snv_annotation import (INVARIANT, INSERTION, DELETION, SNP,
                                         COMPLEX, INDEL, calculate_maf_frequency)
from franklin.seq.seqs import get_seq_name
from franklin.snv.snv_filters import SnvNamer
from franklin.snv.snv_annotation import (_allele_count, _get_group,
                                         calculate_snv_kind)
from franklin.utils.misc_utils import OrderedDict
from franklin.utils.cmd_utils import call

class SnvSequenomWriter(object):
    'It writes the snv in the sequenom format'
    def __init__(self, fhand, length=150, maf=None):
        'It initiates the class'
        self.fhand = fhand
        self.num_features = 0
        self._length = length
        self._maf = maf

    @staticmethod
    def snv_in_illumina(snv):
        "it checks if the snp is only for illumina"
        return not snv.qualifiers['filters']['is_variable'][('samples', ('rp_75_59_uc82',), True)]

    def write(self, sequence, selected_snv_location):
        'It writes a seq with the alternative alleles in one position and Ns in the others.'
        start = selected_snv_location - self._length
        end =  selected_snv_location + self._length + 1
        if start < 0:
            start = 0
        if end > len(sequence):
            end = len(sequence)
        sequence = sequence[start: end]

        selected_snv_location -= start
        maf_threshold = self._maf
        prev_seq_end = 0
        seq_to_print = ''
        for snv in sequence.get_features(kind='snv'):
            # snv start and end [start, end[.
            # Correcting the previous sequence slice
            snv_start = snv.location.start.position - start
            snv_end = snv.location.end.position - start
            # join the previous sequence to the sequence to print
            seq_to_print += str(sequence[prev_seq_end:snv_start].seq)
            prev_seq_end = snv_end

            if snv_start == selected_snv_location:
                #subtituir por allelos
                snv_kind = calculate_snv_kind(snv)
                if snv_kind != SNP:
                    msg = "We don't know how to print anything but SNPs"
                    raise NotImplementedError(msg)
                alleles = '/'.join([a[0] for a in snv.qualifiers['alleles'].keys()])
                to_print = '[{0:s}]'.format(alleles)
            else:
                if maf_threshold is not None:
                    snv_maf = calculate_maf_frequency(snv)
                    write_abundant_allele = True if snv_maf > maf_threshold else False
                else:
                    write_abundant_allele = False
                if write_abundant_allele:
                    # most abundant allele
                    to_print = _get_major_allele(snv)
                else:
                    # Ns
                    snv_kind = calculate_snv_kind(snv)
                    if snv_kind == SNP:
                        to_print = _snp_to_iupac(snv, sequence)
                    elif snv_kind in (DELETION, COMPLEX, INDEL):
                        ref_allele = snv.qualifiers['reference_allele']
                        to_print = ref_allele[0] + 'N' * (len(ref_allele) - 1)
                    else:
                        to_print = 'N'

            seq_to_print += to_print
        else:
            seq_to_print += str(sequence[prev_seq_end:end + 1].seq)

        name = sequence.name + '_' + str(selected_snv_location + 1)
        self.fhand.write('>%s\n%s\n' % (name, seq_to_print))
        self.fhand.flush()


    def write_old(self, sequence, position):
        'It does the real write of the features'
        seq = list(sequence.seq)
        # put N in the snps near of the position
        for snv in sequence.get_features(kind='snv'):
            location = snv.location.start.position
            if location == position:
                mysnv = snv

            elif abs(location - position) < self._length:
                genotype = _snv_to_n(snv, sequence, position, maf=self._maf)
                for index, allele in enumerate(genotype):
                    seq[location + index] = allele

        # Calculate sequence limits
        left_limit  = position - self._length
        rigth_limit = position + self._length + 1
        if left_limit < 0:
            left_limit = 0
        if rigth_limit > len(sequence):
            rigth_limit = len(sequence)

        # sequence as a list, adding
        seq_to_print = list(seq[left_limit: position])
        allele, size = _snv_to_string(mysnv, sequence, position)
        seq_to_print.append(allele)

        seq_to_print.extend(list(seq[position + 1 + size - 1: rigth_limit]))
        seq_to_print = ''.join(seq_to_print)

        name = sequence.name + '_' + str(position + 1)
        self.fhand.write('>%s\n%s\n' % (name, seq_to_print))
        self.fhand.flush()

def _snp_to_iupac(snv, seq):
    '''It converts a snp into its iupac code'''

    iupac_code = {'M': ['A', 'C'],
                  'R': ['A', 'G'],
                  'W': ['A', 'T'],
                  'S': ['C', 'G'],
                  'Y': ['C', 'T'],
                  'K': ['G', 'T'],
                  'V': ['A', 'C', 'G'],
                  'H': ['A', 'C', 'T'],
                  'D': ['A', 'G', 'T'],
                  'B': ['C', 'G', 'T'],
                  'N': ['A', 'C', 'G', 'T'],
                  }
    snp = []
    kind = []
    alleles = snv.qualifiers['alleles'].keys()

    for allele in alleles:
        snp.append(allele[0])
        kind.append(allele[1])

    if 1 in kind or 2 in kind:
        return '[COMPLEX]'
    else:
        sorted_snp = sorted(snp)
        for code in iupac_code:
            if iupac_code[code] == sorted_snp:
                return code

    raise ValueError('Error in getting SNP IUPAC code')

def _get_major_allele(snv):
    'It returns the most frequent allele'
    alleles = snv.qualifiers['alleles']
    major_number_reads = None
    most_freq_allele = None
    for allele in alleles:
        number_reads = _allele_count(allele, alleles)
        if major_number_reads is None or major_number_reads < number_reads:
            major_number_reads = number_reads
            most_freq_allele = allele
    return most_freq_allele[0]

def _snv_to_n(snv, sequence, position, maf=None):
    'It returns the n for each snp'
    genotype = []
    for allele, kind in snv.qualifiers['alleles'].keys():
        if kind == SNP and not genotype:
            snv_maf = calculate_maf_frequency(snv)
            if maf and snv_maf > maf:
                genotype = [_get_major_allele(snv)]
            else:
                snp_iupac = _snp_to_iupac(snv, sequence)
                genotype = [snp_iupac]

        elif kind == DELETION:
            len_del = len(allele)
            genotype.extend(['N'] * (len_del - len(genotype)))
        elif kind == INSERTION:
            geno = sequence[position] + len(allele) * 'N'
            if genotype:
                genotype[0] = geno
            else:
                genotype.append(geno)
    return genotype

def _snv_to_string(snv, sequence, position):
    '''it writes the snv in a format with braquets and with all alleles:
    [A/T], [-/ATG]...'''
    snv_kind = calculate_snv_kind(snv)

    reference_allele = snv.qualifiers['reference_allele']
    alleles = [allele[0] for allele in snv.qualifiers['alleles'].keys()  if allele[1] != INVARIANT]
    alleles = set(alleles)
    assert reference_allele.upper() in ['A', 'T', 'C', 'G']

    if snv_kind == SNP:
        to_print = "%s/%s" % (reference_allele, "/".join(alleles))
        size     = 1
    elif snv_kind == COMPLEX:
        raise RuntimeError('Complex type not implemented')
    elif snv_kind == INDEL:
        if len(alleles) > 1:
            raise RuntimeError('INDEL type not implemented')
    elif snv_kind == INSERTION:
        allele = list(alleles)[0]
        size = 1
        to_print = '-/' + allele
    elif snv_kind == DELETION:
        allele = list(alleles)[0]
        size   = len(allele)

        deleted_alleles = sequence[position:position+size]
        to_print = '%s/-' % ''.join(deleted_alleles)

    to_print = '[%s]' % to_print
    return to_print, size

class SnvIlluminaWriter(object):
    'It writes the snv in the illumina genotyper format'
    def __init__(self, fhand, length=60, write_short=False):
        'It initiates the class'
        self.fhand = fhand
        self.num_features = 0
        self._length = length
        self._write_short = write_short

    def write(self, sequence):
        'It does the real write of the features'
        seq_name = get_seq_name(sequence)

        for snv in sequence.get_features(kind='snv'):
            self.num_features += 1
            location  =  snv.location.start.position
            reference_allele = snv.qualifiers['reference_allele']
            snv_name  =  "%s_%d" % (seq_name, location + 1)
            left_limit  = location - self._length
            rigth_limit = location + self._length + 1
            if self._write_short and left_limit < 0:
                left_limit = 0
            if self._write_short and rigth_limit > len(sequence):
                rigth_limit = len(sequence)


            seq_left   = sequence[left_limit: location]
            seq_rigth  = sequence[location + 1: rigth_limit]
            alleles = [allele[0] for allele in snv.qualifiers['alleles'].keys()]
            alleles = set(alleles)
            alleles.add(reference_allele)
            alleles_str = "[" + "/".join(alleles) + "]"
            illum_str = '%s,SNP,%s\n' % (snv_name,
                                     seq_left.seq + alleles_str + seq_rigth.seq)
            self.fhand.write(illum_str)
        self.fhand.flush()
def compress_and_index_vcf(vcf_fpath):
    '''It indexes the vcf file using tabix and bgzip. the indexes file will be
    vcf_filename.gz
    '''
    cmd = ['bgzip', '-f', vcf_fpath]
    call(cmd, raise_on_error=True)

    cmd = ['tabix', '-p', 'vcf', '-f', '{0:s}.gz'.format(vcf_fpath)]
    call(cmd, raise_on_error=True)

#http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
class VariantCallFormatWriter(object):
    'It writes variant call format files for the snvs.'
    def __init__(self, fhand, reference_name, grouping=None):
        'It inits the class'
        # The fhand is as it arrives
        open(fhand.name, 'w')
        self.fhand = open(fhand.name, 'a')
        self._namer = SnvNamer()

        self._temp_fhand = NamedTemporaryFile(mode='a')
        self._filter_descriptions = {}
        self._header = []
        if grouping is None:
            grouping = 'read_groups'
        self._genotype_grouping_key = grouping
        self._genotype_groups = OrderedDict()
        self._get_pre_header(reference_name)
        self.num_features = 0

    def _get_pre_header(self, reference_name):
        'It writes the header of the vcf file'
        header = self._header
        header.append('##fileformat=VCFv4.1')
        header.append('##fileDate=%s' %
                                      datetime.date.today().strftime('%Y%m%d'))
        header.append('##source=franklin')
        header.append('##reference=%s' % reference_name)
        header.append('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">')
        header.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
        header.append('##INFO=<ID=RC,Number=A,Type=Integer,Description="Read count of the alt alleles">')
        header.append('##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">')
        header.append('##INFO=<ID=BQ,Number=1,Type=Float,Description="RMS Base Quality">')
        header.append('##INFO=<ID=GC,Number=.,Type=String,Description="Genotype Counts: Num. genotypes in which every alleles has been detected">')
        header.append('##INFO=<ID=GP,Number=1,Type=String,Description="Genotype polimorphism">')
        header.append('##INFO=<ID=EZ,Number=1,Type=String,Description="CAP enzymes">')
        header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Read group Genotype">')
        header.append('##FORMAT=<ID=EC,Number=.,Type=Integer,Description="Allele count for the ref and alt alleles in the order listed">')

    def close(self):
        'It merges the header and the snv data'
        # Append the data spec  to the header
        fhand = self.fhand
        self._add_filters_to_header()
        line_items = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                      'INFO', 'FORMAT']
        line_items.extend([group.upper() for group in self._genotype_groups.keys()])
        num_items_per_line = len(line_items)
        self._header.append('%s\n' % '\t'.join(line_items))
        fhand.write('\n'.join(self._header))
        for line in open(self._temp_fhand.name):
            #fix the missing genotype groups
            line = line.strip()
            line = line.split()
            num_items = len(line)
            line.extend(['.:.'] * (num_items_per_line - num_items))
            line = '\t'.join(line)
            #\fix the missing genotype groups
            fhand.write(line + '\n')
        fhand.flush()

    def _add_filters_to_header(self):
        'It adds the used filter tag to the header'
        for name, desc in self._filter_descriptions.values():
            filter_desc = '##FILTER=<ID=%s,Description="%s">' % (name, desc)
            self._header.append(filter_desc)

    def write(self, sequence):
        'It writes the snvs present in the given sequence as SeqFeatures'
        for snv in sequence.get_features(kind='snv'):
            self.num_features += 1
            self._write_snv(sequence, snv)

    @staticmethod
    def _create_alternative_alleles(alleles):
        'It returns the ALT part on the vcf'
        str_alleles = []
        alternative_alleles = []
        for allele in alleles:
            kind = allele[1]
            if kind == INVARIANT:
                continue

            str_allele = allele[0].replace('-', '')
            str_alleles.append(str_allele)
            alternative_alleles.append(allele)
        if str_alleles:
            str_alleles = ','.join(str_alleles)
        else:
            str_alleles = '.'
        return str_alleles, alternative_alleles

    def _create_filters(self, qualifiers):
        'It returns the FILTER part on the vcf'
        filter_strs = []
        if 'filters' not in qualifiers:
            return '.'
        for name, filters_data in qualifiers['filters'].items():
            for parameters, result in filters_data.items():
                if not result:
                    continue
                short_name, description = self._namer.get_filter_description(
                                                      name, parameters,
                                                      self._filter_descriptions)
                short_name = short_name.upper()
                filter_strs.append(short_name)
                self._filter_descriptions[name, parameters] = (short_name,
                                                               description)
        if not filter_strs:
            return 'PASS'
        else:
            return ';'.join(filter_strs)

    def _create_info(self, qualifiers, alternative_alleles):
        'It creates the INFO bit on the vcf'
        toprint_items = []

        alleles = qualifiers['alleles']

        #RC allele count in genotypes, for each ALT allele, in the same order as
        #listed
        acounts = [] #allele_count
        for allele in alternative_alleles:
            acount = _allele_count(allele, alleles, group_kind='read_groups')
            acounts.append(acount)
        if acounts:
            toprint_items.append('RC=%s' % ','.join(map(str, acounts)))

        #AF allele frequency for each ALT allele in the same order as listed:
        reference_allele = qualifiers['reference_allele'], INVARIANT
        if reference_allele in alleles:
            ref_count = _allele_count(reference_allele, alleles,
                                      group_kind='read_groups')
        else:
            ref_count = 0
        total_count = float(sum(acounts) + ref_count)
        afreqs = [acount / total_count for acount in acounts]
        if afreqs:
            toprint_items.append('AF=%s' % ','.join(map(lambda x: '%.1f' % x,
                                                        afreqs)))

        #MQ RMS mapping quality, e.g. MQ=52
        #BQ RMS base quality at this position
        for kind, strfmt in (('mapping_quality', 'MQ=%.2f'),
                             ('quality', 'BQ=%.2f')):
            qual = qualifiers[kind]
            if qual is not None:
                toprint_items.append(strfmt % qual)

        #genotype count
        #we count in how many genotypes every allele has been found.
        allele_counts = self._allele_count_by_group(alleles=alleles,
                                           reference_allele=reference_allele[0],
                                        alternative_alleles=alternative_alleles,
                                          read_groups=qualifiers['read_groups'],
                                          count_reads=False)
        #if some allele is missing there are 0 counts of it
        n_als = max(allele_counts.keys()) + 1
        genotype_counts = [(al, len(allele_counts.get(al, [])))for al in range(n_als)]

        #now we print
        counts = [str(count[1]) for count in genotype_counts]
        toprint_items.append('GC=%s' % (','.join(counts)))
        #genotype polymorphism
        #1 - (number_groups_for_the_allele_with_more_groups / number_groups)
        number_of_groups = sum([count[1] for count in genotype_counts])

        genotype_polymorphism = 1 - genotype_counts[0][1] / float(number_of_groups)
        toprint_items.append('GP=%.2f' % genotype_polymorphism)

        #cap enzymes
        if 'cap_enzymes' in qualifiers and qualifiers['cap_enzymes']:
            to_print = 'EZ=%s' % ','.join(qualifiers['cap_enzymes'])
            toprint_items.append(to_print)

        if toprint_items:
            return ';'.join(toprint_items)
        else:
            return '.'

    @staticmethod
    def _create_quality(alleles, alternative_alleles):
        '''It returns the quality for this snv

        QUAL phred-scaled quality score for the assertion made in ALT. i.e. give
        -10log_10 prob(call in ALT is wrong). If ALT is "." (no variant) then
        this is -10log_10 p(variant), and if ALT is not "." this is -10log_10
        p(no variant). High QUAL scores indicate high confidence calls.
        Although traditionally people use integer phred scores, this field is
        permitted to be a floating point so to enable higher resolution for low
        confidence calls if desired. (Numeric, Missing Value: -1)'''

        if alternative_alleles:
            phreds = [alleles[allele]['quality'] for allele in alternative_alleles]
            if len(phreds) == 1:
                phred = phreds[0]
            else:
                inv_phred = lambda phred: math.pow(10, (-phred / 10))
                probs = map(inv_phred, phreds[:2])
                prob = probs[0] * probs[1]
                phred = -10 * math.log10(prob)
        else:
            phred = alleles.values()[0]['quality']
        return '%i' % phred

    @staticmethod
    def _numbers_for_alleles(reference_allele, alternative_alleles):
        'It returns a key with the numbers for the alleles'
        #a map from alleles to allele index (0 for reference, etc)
        alleles_index = [(reference_allele, INVARIANT)]
        alleles_index.extend(alternative_alleles)
        alleles_index = dict(zip(alleles_index, range(len(alleles_index))))
        return alleles_index

    def _allele_count_by_group(self, alleles, reference_allele,
                               alternative_alleles, read_groups, count_reads):
        '''It returns the allele counts by group

        It can answer to two questions:
            - How many times have the allele been found in every group?
              (count_reads = True)
              It returns a dict indexed by group and the alleles
              (with the vcf number coding)
            - in which groups the allele have been found?
              (count_reads= False)
              It returns a dict indexed by the alleles
              (with the vcf number coding)
        It takes into account the grouping_key (read_group, sample or library)
        It requires the dict with the information about the read_groups.
        '''
        #a map from alleles to allele index (0 for reference, etc)

        alleles_index = self._numbers_for_alleles(reference_allele,
                                                  alternative_alleles)
        grouping_key = self._genotype_grouping_key

        alleles_by_group = {}
        for allele, allele_info in alleles.items():
            #we need the index for the allele
#            if allele not in alleles_index:
#                continue

            try:
                allele_index = alleles_index[allele]
            except KeyError:
                print 'allele', allele
                print "allele index", alleles_index
                print 'ref_allele', reference_allele
                print 'alternative_alleles', alternative_alleles
                print 'alleles', alleles
                raise

            for read_group in allele_info['read_groups']:
                group = _get_group(read_group, grouping_key, read_groups)
                if group not in self._genotype_groups:
                    self._genotype_groups[group] = True
                if count_reads:
                    if group not in alleles_by_group:
                        alleles_by_group[group] = {}
                    if allele_index not in alleles_by_group[group]:
                        alleles_by_group[group][allele_index] = 0
                    count = allele_info['read_groups'][read_group]
                    alleles_by_group[group][allele_index] += count
                else:
                    if allele_index not in alleles_by_group:
                        alleles_by_group[allele_index] = set()
                    alleles_by_group[allele_index].add(group)

        return alleles_by_group

    def _create_genotypes(self, qualifiers, alternative_alleles):
        'It returns the genotype section for this snv'

        alleles = qualifiers['alleles']
        reference_allele = qualifiers['reference_allele']
        read_groups = qualifiers['read_groups']

        items = []
        #the format
        items.append('GT:EC')

        #a map from alleles to allele index (0 for reference, etc)
        alleles_index = self._numbers_for_alleles(reference_allele,
                                                  alternative_alleles)

        #now we need the alleles for every sample

        alleles_by_group = self._allele_count_by_group(alleles=alleles,
                                              reference_allele=reference_allele,
                                        alternative_alleles=alternative_alleles,
                                                        read_groups=read_groups,
                                                               count_reads=True)

        #now we can build the info for every sample
        for group in self._genotype_groups.keys():
            allele_counts = alleles_by_group.get(group, {None: None})
            alleles, counts = [], []
            for allele, count in allele_counts.items():
                allele = str(allele) if allele is not None else '.'
                alleles.append(str(allele))
                count = str(count) if count is not None else '.'
                counts.append(count)
            #print group, alleles, counts
            mix_genotype = '|'.join(alleles)
            allele_counts = ','.join(counts)
            items.append('%s:%s' % (mix_genotype, allele_counts))
        #print 'result', ' '.join(items)

        return '\t'.join(items)

    def _write_snv(self, sequence, snv):
        'Given an snv feature it writes a line in the vcf'
        items = [] #items to write
        items.append(get_seq_name(sequence))
        items.append(str(int(snv.location.start.position) + 1))
        id_ = snv.id
        if id_ == "<unknown id>":
            id_ = '.'
        items.append(id_)
        qualifiers = snv.qualifiers
        ref_seq = qualifiers['reference_allele'].replace('-', '')
        items.append(ref_seq)
        toprint_af, alternative_alleles = self._create_alternative_alleles(
                                                          qualifiers['alleles'])
        items.append(toprint_af)
        items.append(self._create_quality(qualifiers['alleles'],
                                          alternative_alleles))
        filters = self._create_filters(qualifiers)
        items.append(filters)
        try:
            items.append(self._create_info(qualifiers, alternative_alleles))
        except KeyError:
            print 'sequence', get_seq_name(sequence)
            print 'position', str(int(snv.location.start.position))
            raise

        items.append(self._create_genotypes(qualifiers, alternative_alleles))

        self._temp_fhand.write('%s\n' % '\t'.join(items))
        self._temp_fhand.flush()
