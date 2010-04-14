'''
Created on 25/03/2010

@author: peio
'''
import unittest
from StringIO import StringIO
from franklin.seq.seqs import Seq, SeqWithQuality
from franklin.seq.writers import GffWriter, SsrWriter
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

class WriterTest(unittest.TestCase):
    'It test all writers for seqs'

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
        fhand = StringIO()
        ssrwriter = SsrWriter(fhand)
        ssrwriter.write(seq)
        result = fhand.getvalue()
        print result
        assert "seq1\t0\t29\t29\t27\ttrinucleotide\tATC" in result
    @staticmethod
    def test_gff_writer():
        'it tests gff writer'
        srr_feature = SeqFeature(FeatureLocation(ExactPosition(0),
                                                 ExactPosition(29)),
                                                 type='microsatellite',
                                                 qualifiers={'score': 27,
                                                       'type': 'trinucleotide',
                                                       'unit': 'ATC'})
        intron_feature = SeqFeature(FeatureLocation(ExactPosition(478),
                                                    ExactPosition(478)),
                                                    type='intron',
                                          qualifiers={'genomic_db': 'pathtodb'})
        orf_feature = SeqFeature(FeatureLocation(ExactPosition(61),
                                                 ExactPosition(471)),
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

        fhand = StringIO()
        gff_writer = GffWriter(fhand, source='454_roche')
        gff_writer.write(seq)
        gff = fhand.getvalue()
        assert "description=Some desc" in gff
        assert "Ontology_term=GO:0019253" in gff
        assert 'ID=seq1_microsatellite_1'in gff
        assert 'ID=seq1_orf_1;' in gff
        assert 'ID=seq1_intron_1' in gff
        assert 'orthologs=arabidopsis:ara1' in gff
        assert 'name=seq1_snv_1' in gff





if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_gff_writer']
    unittest.main()
