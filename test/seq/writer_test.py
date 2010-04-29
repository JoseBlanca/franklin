'''
Created on 25/03/2010

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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

from Bio.Alphabet import ProteinAlphabet, DNAAlphabet
import unittest
from StringIO import StringIO
from franklin.seq.seqs import Seq, SeqWithQuality
from franklin.seq.writers import GffWriter, SsrWriter, OrfWriter, OrthologWriter
from tempfile import NamedTemporaryFile
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

class WriterTest(unittest.TestCase):
    'It test all writers for seqs'


    @staticmethod
    def test_ortologth_writer():
        'It test ssr writer'
        annotations = {'arabidpsis-orthologs':['ara1', 'ara2'],
                       'melon-orthologs':['melo1', 'melo2']}

        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=[], annotations=annotations)
        fhand = NamedTemporaryFile(mode='a')

        ortologwriter = OrthologWriter(fhand)
        ortologwriter.write(seq)
        seq_result = open(fhand.name).read()
        assert  'ara1,ara2' in seq_result
        assert  'melo1,melo2' in seq_result


    @staticmethod
    def test_orf_writer():
        'It test ssr writer'
        qualifiers = {'pep': Seq('LHPFSHPPXWPLX', ProteinAlphabet()),
                      'dna': Seq('ATGGCTTCATCCATTCTCTCATCCGCCG', DNAAlphabet()),
                      'strand': 'reverse'}
        orf_feature = SeqFeature(FeatureLocation(ExactPosition(61),
                                                 ExactPosition(471)),
                                                 type='orf',
                                                 qualifiers=qualifiers)

        features = [orf_feature]
        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=features, annotations={})
        fhand = NamedTemporaryFile(mode='a')
        pep_fhand = NamedTemporaryFile(mode='a')
        orfwriter = OrfWriter(fhand, pep_fhand)
        orfwriter.write(seq)
        seq_result = open(fhand.name).read()
        pep_result = open(pep_fhand.name).read()
        assert 'ATGGCTTCATCCATTCTCTCATCCGCCG' in seq_result
        assert 'LHPFSHPPXWPLX' in pep_result



    @staticmethod
    def test_ssr_writer():
        'It test ssr writer'
        srr_feature = SeqFeature(FeatureLocation(ExactPosition(0),
                                                 ExactPosition(29)),
                                                 type='microsatellite',
                                                 qualifiers={'score': 27,
                                                       'type': 'trinucleotide',
                                                       'unit': 'ATC'})
        features = [srr_feature]
        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=features, annotations={})
        fhand = NamedTemporaryFile(mode='a')
        ssrwriter = SsrWriter(fhand)
        ssrwriter.write(seq)
        result = open(fhand.name).read()
        assert "seq1\t0\t29\t30\t27\ttrinucleotide\tATC" in result

    @staticmethod
    def test_gff_writer():
        'it tests gff writer'
        srr_feature = SeqFeature(FeatureLocation(ExactPosition(0),
                                                 ExactPosition(29)),
                                                 type='microsatellite',
                                                 qualifiers={'score': 27,
                                                       'type': 'trinucleotide',
                                                       'unit': 'ATC'})
        intron_feature = SeqFeature(FeatureLocation(ExactPosition(30),
                                                    ExactPosition(30)),
                                                    type='intron',
                                          qualifiers={'genomic_db': 'pathtodb'})
        orf_feature = SeqFeature(FeatureLocation(ExactPosition(2),
                                                 ExactPosition(35)),
                                                 type='orf',
                                        qualifiers={'pep': Seq('MASSILSSAXVA'),
                                                    'dna': Seq('ATGGCTTCATCC'),
                                                    'strand':'forward'})
        alleles = {('A', 0): {'read_names':['r1']}}
        snv = SeqFeature(location=FeatureLocation(3, 3), type='snv',
                          qualifiers={'alleles':alleles})
        features = [srr_feature, intron_feature, orf_feature, snv]

        annotations={'GOs': ['GO:0019253', 'GO:0016984', ],
                     'arabidopsis-orthologs':['ara1', 'ara2'],
                     'melo-orthologs':['mel1', 'mel2']}

        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=features,
                             annotations=annotations)

        srr_feature = SeqFeature(FeatureLocation(ExactPosition(0),
                                                 ExactPosition(29)),
                                                 type='microsatellite',
                                                 qualifiers={'score': 27,
                                                       'type': 'trinucleotide',
                                                       'unit': 'ATC'})
        intron_feature = SeqFeature(FeatureLocation(ExactPosition(34),
                                                    ExactPosition(34)),
                                                    type='intron',
                                          qualifiers={'genomic_db': 'pathtodb'})
        orf_feature = SeqFeature(FeatureLocation(ExactPosition(10),
                                                 ExactPosition(15)),
                                                 type='orf',
                                        qualifiers={'pep': Seq('MASSILSSAXVA'),
                                                    'dna': Seq('ATGGCTTCATCC'),
                                                    'strand':'forward'})
        alleles = {('A', 0): {'read_names':['r1']}}
        snv = SeqFeature(location=FeatureLocation(18, 18), type='snv',
                          qualifiers={'alleles':alleles})
        features = [srr_feature, intron_feature, orf_feature, snv]

        annotations={'GOs': ['GO:0019253', 'GO:0016984', ],
                     'arabidopsis-orthologs':['ara1', 'ara2'],
                     'melo-orthologs':['mel1', 'mel2']}

        seq2 = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq2', name='seq2', description='Some desc',
                             dbxrefs=[], features=features,
                             annotations=annotations)
        fhand = StringIO()
        gff_writer = GffWriter(fhand, source='454_roche')
        gff_writer.write(seq)
        gff_writer.write(seq2)
        gff = fhand.getvalue()

        assert "description=Some desc" in gff
        assert "Ontology_term=GO:0019253" in gff
        assert 'ID=seq1_microsatellite_1'in gff
        assert 'ID=seq1_orf_1;' in gff
        assert 'ID=seq1_intron_1' in gff
        assert 'orthologs=arabidopsis:ara1' in gff
        assert 'name=seq1_snv_1' in gff


        orf_feature = SeqFeature(FeatureLocation(ExactPosition(3),
                                                 ExactPosition(5)),
                                                 type='orf',
                                        qualifiers={'pep': Seq('MASSILSSAXVA'),
                                                    'dna': Seq('ATGGCTTCATCC'),
                                                    'strand':'forward'})

        seq1 = SeqWithQuality(seq=Seq('CTTCATCCAT'),
                             id='seq1', name='seq1', description='Some desc',
                             dbxrefs=[], features=[orf_feature],
                             annotations=annotations)
        seq2 = SeqWithQuality(seq=Seq('CTTCATCCAT'),
                             id='seq2', name='seq2', description='Some desc',
                             dbxrefs=[], features=[orf_feature],
                             annotations=annotations)
        fhand = StringIO()
        gff_writer = GffWriter(fhand, source='454_roche')
        gff_writer.write(seq1)
        gff_writer.write(seq2)
        gff = fhand.getvalue()
        print gff


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_gff_writer']
    unittest.main()
