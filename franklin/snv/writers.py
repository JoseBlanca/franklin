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

from franklin.snv.snv_annotation import INVARIANT, INSERTION, DELETION, SNP
from franklin.seq.seqs import get_seq_name
from franklin.snv.snv_filters import get_filter_description
from franklin.utils.misc_utils import OrderedDict

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

#http://1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf3.3
class VariantCallFormatWriter(object):
    'It writes variant call format files for the snvs.'
    def __init__(self, fhand, reference_name, grouping=None):
        'It inits the class'
        # The fhand is as it arrives
        open(fhand.name, 'w')
        self.fhand = open(fhand.name, 'a')

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
        header.append('##format=VCFv3.3')
        header.append('##fileDate=%s' %
                                      datetime.date.today().strftime('%Y%m%d'))
        header.append('##source=franklin')
        header.append('##reference=%s' % reference_name)
        header.append('##INFO=NS,1,Integer,"Number of Samples With Data"')
        header.append('##INFO=AF,-1,Float,"Allele Frequency"')
        header.append('##INFO=AC,-1,Integer,"Allele Count"')
        header.append('##INFO=MQ,-1,Float,"RMS Mapping Quality"')
        header.append('##INFO=BQ,-1,Float,"RMS Base Quality"')
        header.append('##INFO=GC,.,String,"Genotype counts for alleles"')
        header.append('##FORMAT=RG,.,String,"Read group Genotype genotype"')
        header.append('##FORMAT=AC,.,String,"Allele count"')

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
            filter_desc = '##FILTER=%s,"%s"' % (name, desc)
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
            if kind == SNP:
                str_allele = allele[0]
            elif kind == INSERTION:
                str_allele = 'I%s' % allele[0]
            elif kind == DELETION:
                str_allele = 'D%d' % len(allele[0])
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
                short_name, description = get_filter_description(name,
                                                                 parameters,
                                                      self._filter_descriptions)
                short_name = short_name.upper()
                filter_strs.append(short_name)
                self._filter_descriptions[name, parameters] = (short_name,
                                                               description)
        if not filter_strs:
            return '.'
        else:
            return ';'.join(filter_strs)

    def _create_info(self, qualifiers, alternative_alleles):
        'It creates the INFO bit on the vcf'
        toprint_items = []

        alleles = qualifiers['alleles']

        allele_count = lambda al: sum(alleles[al]['read_groups'].values())

        #AC allele count in genotypes, for each ALT allele, in the same order as
        #listed
        acounts = [] #allele_count
        for allele in alternative_alleles:
            acounts.append(allele_count(allele))
        if acounts:
            toprint_items.append('AC=%s' % ','.join(map(str, acounts)))

        #AF allele frequency for each ALT allele in the same order as listed:
        reference_allele = qualifiers['reference_allele'], INVARIANT
        if reference_allele in alleles:
            ref_count = allele_count(reference_allele)
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
        grouping_key = self._genotype_grouping_key
        allele_numb_coding = self._numbers_for_alleles(reference_allele[0],
                                                       alternative_alleles)
        #we count the number of times every allele appears in a group (sample)
        genotype_counts = []
        for allele, allele_info in qualifiers['alleles'].items():
            allele = allele_numb_coding[allele]
            genotype_counts.append((allele,
                                    len(set(allele_info[grouping_key]))))
        #we have to sort by the number of groups
        genotype_counts = sorted(genotype_counts, lambda x, y: y[1] - x[1])
        #now we print
        alleles = [str(count[0]) for count in genotype_counts]
        counts = [str(count[1]) for count in genotype_counts]
        toprint_items.append('GC=%s:%s' % ('|'.join(alleles), ','.join(counts)))

        #genotype polymorphism
        #1 - (number_groups_for_the_allele_with_more_groups) / number_groups
        number_of_groups = sum([count[1] for count in genotype_counts])
        genotype_polymorphism = 1 - genotype_counts[0][1] / float(number_of_groups)
        toprint_items.append('GP=%.2f' % genotype_polymorphism)

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

    def _create_genotypes(self, alleles, reference_allele, alternative_alleles):
        'It returns the genotype section for this snv'

        items = []
        #the format
        items.append('RG:AC')
        #a map from alleles to allele index (0 for reference, etc)
        alleles_index = self._numbers_for_alleles(reference_allele,
                                                  alternative_alleles)

        #now we need the alleles for every sample
        grouping_key = self._genotype_grouping_key
        alleles_by_group = {}
        for allele, allele_info in alleles.items():
            #we need the index for the allele
            allele_index = alleles_index[allele]
            for group in allele_info[grouping_key]:
                if group not in self._genotype_groups:
                    self._genotype_groups[group] = True
                if group not in alleles_by_group:
                    alleles_by_group[group] = {}
                if allele_index not in alleles_by_group[group]:
                    alleles_by_group[group][allele_index] = 0
                alleles_by_group[group][allele_index] += 1

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
        ref_seq = qualifiers['reference_allele']
        items.append(ref_seq)
        toprint_af, alternative_alleles = self._create_alternative_alleles(
                                                          qualifiers['alleles'])
        items.append(toprint_af)
        items.append(self._create_quality(qualifiers['alleles'],
                                          alternative_alleles))
        filters = self._create_filters(qualifiers)
        items.append(filters)
        items.append(self._create_info(qualifiers, alternative_alleles))
        items.append(self._create_genotypes(qualifiers['alleles'],
                                            qualifiers['reference_allele'],
                                            alternative_alleles))
        self._temp_fhand.write('%s\n' % '\t'.join(items))
        self._temp_fhand.flush()
