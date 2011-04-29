'''
Created on 26/10/2009

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

import unittest, os
from StringIO import StringIO
from tempfile import NamedTemporaryFile

from franklin.gff import (_add_dbxref_to_feature, GffFile, write_gff,
                          METADATA, FEATURE, create_dbxref_feature_mapper,
                          SeqGffWriter, create_go_annot_mapper,
                          create_feature_type_filter)
from franklin.utils.misc_utils import TEST_DATA_DIR
from franklin.seq.seqs import Seq, SeqWithQuality
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition

GFF_CONTENT = '''##gff-version 3
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tID=gene00001;Name=EDEN
##FASTA\n>hola\nGATA\n'''

class GffTest(unittest.TestCase):
    'It tests the GFF class'

    def _create_gff_file(self):
        'It creates a temporary gff for testing, it returns the fhand'

        gff_fhand = NamedTemporaryFile()
        gff_fhand.write(GFF_CONTENT)
        gff_fhand.flush()
        return gff_fhand

    def test_read_gff(self):
        'It opens a GFF for reading'
        fhand = self._create_gff_file()
        gff = GffFile(fpath=fhand.name)
        assert gff.version == '3'
        items = list(gff.items)
        assert items[0][0] == METADATA
        assert items[0][1] == 'sequence-region ctg123 1 1497228'
        assert items[1][0] == FEATURE
        assert items[1][1] == {'end': 9000, 'name': 'EDEN', 'start': 1000,
                               'source': '.', 'seqid': 'ctg123', 'phase': '.',
                               'attributes': {'ID': 'gene00001',
                                              'Name': 'EDEN'},
                               'score': '.', 'type': 'gene', 'id': 'gene00001',
                               'strand': '.'}
        seqs = list(items[2][1])
        assert seqs[0].seq == 'GATA'

    def test_write_gff(self):
        'It writes a gff file'
        gff_out_fhand = NamedTemporaryFile()

        gff_out = GffFile(fpath = gff_out_fhand.name, mode='w')
        gff_in_fhand = self._create_gff_file()
        gff = GffFile(fpath=gff_in_fhand.name)
        for kind, item in gff.items:
            gff_out.write(kind, item)
        gff_out.flush()

        result = open(gff_out_fhand.name).read()
        assert result == GFF_CONTENT

    @staticmethod
    def test_items_in_gff():
        'It gets the items in a gff file'
        gff = GffFile(os.path.join(TEST_DATA_DIR, 'map_fis.gff3'))
        features = [item[1] for item in gff.items if item[0] == FEATURE]
        assert len(features) == 99
        assert features[1]['name'] == 'ctg,0'
        assert features[1]['source'] == 'F=PC'
        assert features[1]['attributes']['Name'] == 'ctg,0'
        assert features[1]['id'] == 'ctg 0'
        assert features[98]['name'] == 'Cm13_B04'

    def test_gff_mappers(self):
        'It test that the mappers in gff are runs'
        database = 'database'
        relations = 'gene00001\tacc1\n'
        rels_fhand = StringIO(relations)
        dbxref_mapper = create_dbxref_feature_mapper(database, rels_fhand)

        fhand = self._create_gff_file()
        gff   = GffFile(fpath=fhand.name, feature_mappers=[dbxref_mapper])
        items = list(gff.items)
        assert items[1][1]['attributes']['Dbxref'] == 'database:acc1'

    def test_gff_filters(self):
        'It test that we can use filters for features in gff'
        types = ['gene']
        fhand = self._create_gff_file()
        type_filter = create_feature_type_filter(types)
        gff   = GffFile(fpath=fhand.name, feature_filters=[type_filter])
        items = list(gff.items)
        assert items[1][1]['id'] == 'gene00001'



        types = ['contig']
        fhand = self._create_gff_file()
        type_filter = create_feature_type_filter(types)
        gff   = GffFile(fpath=fhand.name, feature_filters=[type_filter])
        items = list(gff.items)
        assert len(items) == 2

class GffMappersTest(unittest.TestCase):
    'It test the mappers in GffFile'

    @staticmethod
    def test_dbxref_feature_mapper():
        'It tests the dbxref feature mapper'

        database = 'database'
        relations = 'ctg0\tacc1\n'
        rels_fhand = StringIO(relations)
        feature = {'end': 140722177, 'name': 'ctg0', 'start': 1,
                   'source': 'F=PC', 'seqid': 'Chrctg0', 'phase': '.',
                   'attributes': {'ID': 'ctg0', 'Name': 'ctg,0'},
                   'score': '.', 'type': 'contig', 'id':'ctg0', 'strand': '.'}
        mapper = create_dbxref_feature_mapper(database, rels_fhand)
        mapper(feature)

        # add an already given database
        mapper(feature)
        assert feature['attributes']['Dbxref'] == 'database:acc1'

        # add a second dbxref
        mapper = create_dbxref_feature_mapper('database2', rels_fhand)
        mapper(feature)
        assert feature['attributes']['Dbxref'] == 'database2:acc1,database:acc1'

    @staticmethod
    def test_go_annot_mapper():
        'it test s the go_annot_mapepr'
        go_annot ='''MELO3A000001P1\tGO:0016023\tprotein gi
MELO3A000001P1\tGO:0006950\tprotein gi'''
        annot_fhand = StringIO(go_annot)
        feature = {'end': 140722177, 'name': 'MELO3A000001P1', 'start': 1,
                   'source': 'F=PC', 'seqid': 'Chrctg0', 'phase': '.',
                   'attributes': {'ID': 'MELO3A000001P1', 'Name': 'MELO3A000001P1'},
                   'score': '.', 'type': 'contig', 'id':'MELO3A000001P1',
                   'strand': '.'}
        mapper = create_go_annot_mapper(annot_fhand)
        mapper(feature)
        assert feature['attributes']['Ontology_term'] == 'GO:0016023,GO:0006950'

        #wit already go terms
        feature = {'end': 140722177, 'name': 'MELO3A000001P1', 'start': 1,
                   'source': 'F=PC', 'seqid': 'Chrctg0', 'phase': '.',
                   'attributes': {'ID': 'MELO3A000001P1', 'Name': 'MELO3A000001P1',
                                  'Ontology_term':'GO:0016023'},
                   'score': '.', 'type': 'contig', 'id':'MELO3A000001P1',
                   'strand': '.'}
        mapper(feature)
        assert feature['attributes']['Ontology_term'] == 'GO:0016023,GO:0006950'

        #feature without gos
        feature = {'end': 140722177, 'name': 'MELO3A000001P2', 'start': 1,
                   'source': 'F=PC', 'seqid': 'Chrctg0', 'phase': '.',
                   'attributes': {'ID': 'MELO3A000001P2',
                                  'Name': 'MELO3A000001P2'},
                   'score': '.', 'type': 'contig', 'id':'MELO3A000001P2',
                   'strand': '.'}
        mapper(feature)
        assert 'Ontology_term' not in feature['attributes']

class GffFilterTest(unittest.TestCase):
    'It test the mappers in GffFile'

    @staticmethod
    def test_feature_type_filter():
        'It tests the dbxref feature mapper'

        feature = {'end': 140722177, 'name': 'ctg0', 'start': 1,
                   'source': 'F=PC', 'seqid': 'Chrctg0', 'phase': '.',
                   'attributes': {'ID': 'ctg0', 'Name': 'ctg,0'},
                   'score': '.', 'type': 'contig', 'id':'ctg0', 'strand': '.'}

        type_filter = create_feature_type_filter(['contig'])
        assert filter(type_filter, [feature])

        type_filter = create_feature_type_filter(['gene'])
        assert not filter(type_filter, [feature])

class GffOutTest(unittest.TestCase):
    'It tests the gff writer'

    @staticmethod
    def test_simple_output():
        'We can write a simple gff3 file'
        feat1 = {'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 #'id':    'gene00001',
                 #'name':  'EDEN',
                 'attributes':{'ID':'gene00001', 'Name':'EDEN'}
                 }
        feats = [(METADATA, 'sequence-region ctg123 1 1497228'),
                 (FEATURE, feat1)]
        result = '''##gff-version 3
##sequence-region ctg123 1 1497228
ctg123\t.\tgene\t1000\t9000\t.\t.\t.\tID=gene00001;Name=EDEN\n'''
        outh = NamedTemporaryFile()
        write_gff(outh.name, feats)
        assert result == outh.read()

        feat1 = {'id':'23',
                 'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 'name': 'hola',
                 'attributes' : {'Parent': ['p1', 'p2']}}
        feats = [(FEATURE, feat1)]
        outh = NamedTemporaryFile()
        write_gff(outh.name, feats)
        result = outh.read()
        expected = '##gff-version 3\nctg123\t.\tgene\t1000\t9000\t.\t.\t.\t'
        assert expected in result
        assert 'Name=hola' in result

        #escaping some caracteres
        feat1 = {'id':'23',
                 'seqid': 'ctg123',
                 'type':  'gene',
                 'start': 1000,
                 'end':   9000,
                 'name': 'hola',
                 'attributes':{'Dbxref':'peoi%25l a%20k%s'}}
        feats = [(FEATURE, feat1)]
        outh = NamedTemporaryFile()
        write_gff(outh.name, feats)
        result = outh.read()
        assert 'Dbxref=peoi%25l%20a%20k%25s' in  result

    @staticmethod
    def test_from_2_to_3():
        'It tests that we can go from gff2 to gff3'
        GFF2 = '''Chrctg0\tassembly\tChromosome\t1\t140722177\t.\t.\t.\tSequence "Chrctg0"; Name "Chrctg0"
Chrctg0\tFPC\tcontig\t1\t140722177\t.\t.\t.\tcontig "ctg0"; Name "ctg0"
Chrctg0\tFPC\tBAC\t109076481\t109461505\t.\t.\t.\tBAC "Cm45_J09"; Name "Cm45_J09"; Contig_hit "0"
Chrctg0\tFPC\tBAC\t97189889\t97329153\t.\t.\t.\tBAC "Cm40_O16 3"; Name "Cm40_O16"; Contig_hit "0"
Chrctg0\tFPC\tBAC\t57982977\t58302465\t.\t.\t.\tBAC "Cm22_F20"; Name "Cm22_F20"; Contig_hit "0"
Chrctg0\tFPC\tBAC\t57982978\t58302466\t.\t.\t.\tBAC "Cm22_F20"; Name "Cm22_F20"; Contig_hit "0"
'''
        inh = NamedTemporaryFile()
        inh.write(GFF2)
        inh.flush()
        in_gff = GffFile(inh.name)

        outh = NamedTemporaryFile()
        write_gff(outh.name, in_gff.items)

        result = outh.read()
        assert 'ID=Cm22_F20_2' in result
        assert 'BAC=Cm40_O16%203' in result

    @staticmethod
    def test_add_dbxref():
        'It adds dbxrefs to the features'
        feature = {'seqid': 'ctg123',
                   'type':  'gene',
                   'start': 1000,
                   'end':   9000,
                   'attributes':{'ID':'gene00001', 'Name':'EDEN'}
                   }

        dbxref_db = 'test'
        dbxref_id = 'id100'
        _add_dbxref_to_feature(feature, dbxref_db, dbxref_id)
        assert feature['attributes']['Dbxref'] == 'test:id100'

        feature = {'seqid': 'ctg123',
                   'type':  'gene',
                   'start': 1000,
                   'end':   9000,
                   'attributes':{'ID':'gene00001', 'Name':'EDEN',
                                 'Dbxref':'test2:id101'}
                   }
        _add_dbxref_to_feature(feature, dbxref_db, dbxref_id)
        assert 'test:id100' in feature['attributes']['Dbxref']
        assert 'test2:id101' in feature['attributes']['Dbxref']

class SequenceWriter(unittest.TestCase):
    'It tests that we can write a stream of SeqWithQuality into a file'
    @staticmethod
    def test_gff_writer():
        'it tests the seq to gff writer'
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

        annotations = {'GOs': ['GO:0019253', 'GO:0016984', ],
                     'arabidopsis-orthologs':['ara1', 'ara2'],
                     'melo-orthologs':['mel1', 'mel2']}

        seq = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq1;', name='seq1',
                             description='equal 96%',
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

        annotations = {'GOs': ['GO:0019253', 'GO:0016984', ],
                     'arabidopsis-orthologs':['ara1', 'ara2'],
                     'melo-orthologs':['mel1', 'mel2']}

        seq2 = SeqWithQuality(seq=Seq('CTTCATCCATTCTCTCATCCGCCGNTGTGGCCTTTGN'),
                             id='seq2', name='seq2',
                             dbxrefs=[], features=features,
                             annotations=annotations)
        fhand = NamedTemporaryFile()
        gff_writer = SeqGffWriter(fhand, source='454_roche')
        gff_writer.write(seq)
        gff_writer.write(seq2)
        gff = open(fhand.name).read()

        assert "description=equal%2096%25;" in gff
        assert "Ontology_term=GO:0019253" in gff
        assert 'seq1%3B_microsatellite'in gff
        assert 'ID=seq1%3B_ORF' in gff
        assert 'ID=seq1%3B_intron' in gff
        assert 'orthologs=arabidopsis:ara1' in gff
        assert 'ID=seq1%3B_SNV' in gff

        orf_feature = SeqFeature(FeatureLocation(ExactPosition(3),
                                                 ExactPosition(5)),
                                                 type='orf',
                                        qualifiers={'pep': Seq('MASSILSSAXVA'),
                                                    'dna': Seq('ATGGCTTCATCC'),
                                                    'strand':'forward'})

        seq1 = SeqWithQuality(seq=Seq('CTTCATCCAT'),
                             id='seq1', name='seq1',
                             dbxrefs=[], features=[orf_feature], description='',
                             annotations=annotations)
        seq2 = SeqWithQuality(seq=Seq('CTTCATCCAT'),
                             id='seq2', name='seq2',
                             dbxrefs=[], features=[orf_feature],
                             annotations=annotations)
        fhand = NamedTemporaryFile()
        gff_writer = SeqGffWriter(fhand, source='454_roche')
        gff_writer.write(seq1)
        gff_writer.write(seq2)
        gff = open(fhand.name).read()
        assert 'description' not in gff

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'GffTest.test_gff_filters']
    unittest.main()
