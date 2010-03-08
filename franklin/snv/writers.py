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

import datetime, math
from franklin.snv.snv_annotation import INVARIANT, INSERTION, DELETION, SNP
from franklin.seq.seqs import get_seq_name
from franklin.snv.snv_filters import get_filter_description

#http://1000genomes.org/wiki/doku.php?id=1000_genomes:analysis:vcf3.3

class VariantCallFormatWriter(object):
    'It writes variant call format files for the snvs.'
    def __init__(self, fhand, reference_name):
        'It inits the class'
        open(fhand.name, 'w')
        self._fhand = open(fhand.name, 'a')
        # The fhand is as it arrives
        self._header = []
        self._write_header(reference_name)

    def _write_header(self, reference_name):
        'It writes the header of the vcf file'
        fhand = self._fhand
        self._header.append('##format=VCFv3.3')
        self._header.append('')
        fhand.write('##format=VCFv3.3\n')
        fhand.write('##fileDate=%s\n' %
                                       datetime.date.today().strftime('%Y%m%d'))
        fhand.write('##source=franklin\n')
        fhand.write('##reference=%s\n' % reference_name)

        fhand.write('##INFO=NS,1,Integer,"Number of Samples With Data"\n')
        fhand.write('##INFO=AF,-1,Float,"Allele Frequency"\n')
        fhand.write('##INFO=AC,-1,Integer,"Allele Count"\n')
        fhand.write('##INFO=MQ,-1,Float,"RMS Mapping Quality"\n')
        fhand.write('##INFO=BQ,-1,Float,"RMS Base Quality"\n')
        ##FILTER=q10,"Quality below 10"
        ##FILTER=s50,"Less than 50% of samples have data"

        fhand.write('%s\n' % '\t'.join(('#CHROM', 'POS', 'ID', 'REF', 'ALT',
                                        'QUAL', 'FILTER',  'INFO')))

    def write(self, sequence):
        'It writes the snvs present in the given sequence as SeqFeatures'
        filter_descriptions = {}
        for snv in sequence.get_features(kind='snv'):
            self._write_snv(sequence, snv, filter_descriptions)
        self._fhand.flush()
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

    def _create_filters(self, qualifiers, filter_descriptions):
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
                                                            filter_descriptions)
                filter_strs.append(short_name)
                filter_descriptions[name, parameters] = short_name, description
        if not filter_strs:
            return '.'
        else:
            return ';'.join(filter_strs)

    @staticmethod
    def _root_mean_square(numbers):
        'It returns the root mean square for the given numbers'
        power2 = lambda x: math.pow(x, 2)
        return math.sqrt(sum(map(power2, numbers)) / len(numbers))

    def _create_info(self, qualifiers, alternative_alleles):
        'It creates the INFO bit on the vcf'
        toprint_items = []

        alleles = qualifiers['alleles']

        #AC allele count in genotypes, for each ALT allele, in the same order as
        #listed
        acounts = [] #allele_count
        for allele in alternative_alleles:
            acounts.append(len(alleles[allele]['read_names']))
        if acounts:
            toprint_items.append('AC=%s' % ','.join(map(str, acounts)))

        #AF allele frequency for each ALT allele in the same order as listed:
        reference_allele = qualifiers['reference_allele'], INVARIANT
        if reference_allele in alleles:
            ref_count = len(alleles[reference_allele]['read_names'])
        else:
            ref_count = 0
        total_count = sum(acounts) + ref_count
        afreqs  = [acount/total_count for acount in acounts]
        if afreqs:
            toprint_items.append('AF=%s' % ','.join(map(str, acounts)))

        #MQ RMS mapping quality, e.g. MQ=52
        #BQ RMS base quality at this position
        for kind, strfmt in (('mapping_qualities', 'MQ=%f'),
                             ('qualities', 'BQ=%f')):
            quals = []
            for allele_info in alleles.values():
                quals.extend(allele_info[kind])
            if quals:
                toprint_items.append(strfmt % self._root_mean_square(quals))
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

        phreds = [ alleles[allele]['quality'] for allele in alternative_alleles]
        if len(phreds) == 1:
            phred = phreds[0]
        else:
            inv_phred = lambda phred: math.pow(10, (-phred/10))
            probs = map(inv_phred, phreds[:2])
            prob = probs[0] * probs[1]
            phred = -10 * math.log10(prob)
        return '%i' % phred


    def _write_snv(self, sequence, snv, filter_descriptions):
        'Given an snv feature it writes a line in the vcf'
        items = [] #items to write
        items.append(get_seq_name(sequence))
        items.append(str(snv.location.start.position))
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
        filters = self._create_filters(qualifiers, filter_descriptions)
        #self._add_filter_definition(filters)
        items.append(filters)
        items.append(self._create_info(qualifiers, alternative_alleles))
        self._fhand.write('%s\n' % '\t'.join(items))
