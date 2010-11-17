'''
Created on 29/03/2010

@author: jose
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

import unittest, os, shutil
from os.path import join, exists

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis
from franklin.backbone.analysis import BACKBONE_DIRECTORIES, BACKBONE_BASENAMES
from franklin.utils.cmd_utils import BLAST_TOOL
THREADS = 2

class TestBackboneMapping(unittest.TestCase):
    'It tests the backbone'

    @staticmethod
    def test_mapping_analysis():
        'We can map the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        if BLAST_TOOL == 'blast':
            blastdb_seq = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        else:
            blastdb_seq = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes+')
        snv_filters = {'filter1':{'name':'uniq_contiguous', 'use':True,
                                  'genomic_db':blastdb_seq,
                                  'genomic_seqs_fpath':blastdb_seq},

                       'filter7':{'name':'by_kind', 'use':True,
                                  'kind':'SNP'},
                       'filter12':{'name':'ref_not_in_list', 'use':True,
                                'list_path':os.path.join(DATA_DIR, 'cos_list')},
                       'filter10':{'unique_name': 'variable_in_sm',
                                   'name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['hola']},
                       'filter11':{'unique_name': 'variable_in_adios',
                                   'name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['adios']},
                       'filter13':{'unique_name': 'variable_in_caracola',
                                   'name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['caracola']}, }

        configuration = {'Snvs':{'min_quality':20},
                         'Sam_processing':{'add_default_qualities':True},
                         'snv_filters':snv_filters,
                         'General_settings':{'threads':THREADS}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        solexa = '@seq1\n'
        solexa += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq2\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq14\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq15\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq12\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq13\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa = '@seq16\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'

        sanger = '>seq3\n'
        sanger += 'GATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATAGGAGATACTGA'
        sanger += 'GATTCTGGAATCTCTGAGTTTCTGGGTTCAAGTTGCACTGACCATTGTTGGATTTGTAGA'
        sanger += 'TTGTTTCTTCTTTCATTAGGCATTGATTATGGGTAAATGCGTGGGTACATATAATATATA'
        sanger += 'TCTGTTGAATGCAATTTACACATTGACTGAGGAACAACATGAACATGGCAGCTTTCTCAA'
        sanger += 'AATTGAACCACAGAAGGCTTAAAAGCAAAGTCTTTGGAGAATCAGACTAAGCTTGAGA\n'
        sanger += '>seq4\n'
        sanger += 'TCCATACTTTACTCTATCTCTTTCTGTGTTTGGTAACACGCGAGGATTGGATGATAT'
        sanger += 'CTATATATTCTAATGTGGACTAAAAATGTGTGTGTGTGTATGAAGATGGGAAGCCGGAAG'
        sanger += 'TCATCAAGGAGAAAAGAGAGATAGACGACGAGAAGATGGCGTTGACGTTCAGAGGACTAG'
        sanger += 'AGGGTCATGTGATGGAGAAGTACAAGAAGTATGAGGTTATCTTACAGTTCATTCCCAAGT'
        sanger += 'CGAACGAAGGCTGCGTCTGCAAAGTCACTCTGATATGGGAGAATCGCAACGAAGACTCCC'
        sanger += '>seq5\n'
        sanger += 'TCCATACTTTACTCTATCTCTTTCTGTGTTTGGTAACACGCGAGGATTGGATGATAT'
        sanger += 'CTATATATTCTAAAGTGGACTAAAAATGTGTGTGTGTGTATGAAGATGGGAAGCCGGAAG'
        sanger += 'TCATCAAGGAGAAAAGAGAGATAGACGACGAGAAGATGGCGTTGACGTTCAGAGGACTAG'
        sanger += 'AGGGTCATGTGATGGAGAAGTACAAGAAGTATGAGGTTATCTTACAGTTCATTCCCAAGT'
        sanger += 'CGAACGAAGGCTGCGTCTGCAAAGTCACTCTGATATGGGAGAATCGCAACGAAGACTCCC'

        fpath_sanger = join(clean_reads_dir, 'lb_hola1.pl_sanger.sm_hola.fasta')
        fpath_solexa = join(clean_reads_dir,
                                    'lb_hola2.pl_illumina.sm_hola.sfastq')
        open(fpath_sanger, 'w').write(sanger)
        open(fpath_solexa, 'w').write(solexa)

        fpath_sanger2 = join(clean_reads_dir, 'lb_adios.pl_sanger.fasta')
        fpath_solexa2 = join(clean_reads_dir,
                                    'lb_adios.pl_illumina.sfastq')
        open(fpath_sanger2, 'w').write(sanger)
        open(fpath_solexa2, 'w').write(solexa)

        #the reference
        reference_dir = join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(join(DATA_DIR, 'blast/arabidopsis_genes')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping', silent=True)
        mapping_dir = join(project_dir, 'mapping')
        singular_mapping_dir = sorted(os.listdir(mapping_dir))[0]
        singular_mapping_dir = join(mapping_dir, singular_mapping_dir)
        assert exists(join(singular_mapping_dir, 'bams',
                            'by_readgroup', 'lb_hola2.pl_illumina.sm_hola.bam'))
        result_dir = join(mapping_dir, 'bams')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)

        do_analysis(project_settings=settings_path, kind='merge_bams',
                    silent=True)
        assert exists(join(result_dir, 'merged.0.bam'))
        assert exists(join(result_dir, 'merged.0.bam.bai'))

        #we realign the mapping using GATK
        do_analysis(project_settings=settings_path, kind='realign_bam',
                    silent=True)
        assert exists(join(result_dir, 'merged.1.bam'))

        do_analysis(project_settings=settings_path, kind='mapping_stats',
                    silent=True)
        stats_fname = join(mapping_dir,
                           BACKBONE_DIRECTORIES['mapping_stats'][1],
                           BACKBONE_BASENAMES['statistics_file'])
        result = open(stats_fname).read()
        assert 'Statistics for Coverage for platform sanger' in result

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)
        os.symlink(reference_fpath, join(annot_input_dir, 'reference.fasta'))
        do_analysis(project_settings=settings_path, kind='annotate_snvs',
                    silent=True)
        json_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'reference.0.pickle')
        assert 'snv' in  open(json_fpath).read()

        do_analysis(project_settings=settings_path, kind='filter_snvs',
                    silent=True)
        json_fpath = join(project_dir, BACKBONE_DIRECTORIES['annotation_dbs'],
                          'reference.1.pickle')
        result = open(json_fpath).read()
        #print result
        assert 'snv' in result
        assert 'adios_sanger' in result

        do_analysis(project_settings=settings_path, kind='write_annotations',
                    silent=True)
        vcf_fpath = join(project_dir, 'annotations', 'features',
                         'reference.vcf')
        vcf = open(vcf_fpath).read()
        #print vcf
        assert 'VKS' in vcf
        assert 'VLB1' in vcf
        assert 'VLB2' in vcf
        assert 'VLB3' in vcf
        assert 'AT5G19860.1' in vcf

        do_analysis(project_settings=settings_path, kind='mapping_stats',
                    silent=True)
        stats_dir = join(project_dir, 'mapping', 'bams', 'stats')
        assert exists(join(stats_dir, 'backbone.coverage_illumina.dat'))

        stats_fpath = join(stats_dir, BACKBONE_BASENAMES['statistics_file'])
        result = open(stats_fpath).read()
        expected = 'average: 0.45\nvariance: 1.30\ntotal sequence length: 3941'
        assert expected in result

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_mapping_color():
        'It test the mapping of the mapper with color space'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        if BLAST_TOOL == 'blast':
            blastdb_seq = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        else:
            blastdb_seq = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes+')

        snv_filters = {'filter1':{'name':'uniq_contiguous', 'use':True,
                                  'genomic_db':blastdb_seq,
                                  'genomic_seqs_fpath':blastdb_seq},

                       'filter7':{'name':'by_kind', 'use':True,
                                  'kind':'SNP'},
                       'filter12':{'name':'ref_not_in_list', 'use':True,
                                'list_path':os.path.join(DATA_DIR, 'cos_list')},
                       'filter10':{'unique_name': 'variable_in_sm',
                                   'name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['hola']},
                       'filter11':{'unique_name': 'variable_in_adios',
                                   'name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['adios']},
                       'filter13':{'unique_name': 'variable_in_caracola',
                                   'name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['caracola']}, }

        configuration = {'Snvs':{'min_quality':20},
                         'Sam_processing':{'add_default_qualities':True},
                         'snv_filters':snv_filters,
                         'General_settings':{'threads':THREADS}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)
        shutil.copy(os.path.join(DATA_DIR, 'solid.fastq'),
               os.path.join(clean_reads_dir, 'pl_solid.lb_hola.sm_hola.sfastq'))

        #the reference
        reference_dir = join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(join(DATA_DIR, 'samtools_color/reference')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping', silent=True)
        mapping_dir = join(project_dir, 'mapping')
        singular_mapping_dir = sorted(os.listdir(mapping_dir))[0]
        singular_mapping_dir = join(mapping_dir, singular_mapping_dir)
        assert exists(join(singular_mapping_dir, 'bams',
                            'by_readgroup', 'pl_solid.lb_hola.sm_hola.bam'))
        result_dir = join(mapping_dir, 'bams')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)

        do_analysis(project_settings=settings_path, kind='merge_bams',
                    silent=True)
        assert exists(join(result_dir, 'merged.0.bam'))
        assert exists(join(result_dir, 'merged.0.bam.bai'))

        #we realign the mapping using GATK
        do_analysis(project_settings=settings_path, kind='realign_bam',
                    silent=True)
        assert exists(join(result_dir, 'merged.1.bam'))

        os.chdir('/tmp')
        test_dir.close()

    @staticmethod
    def test_mapping_analysis_boina():
        'We can map the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'

        configuration = {'Sam_processing':{'add_default_qualities':True},
                         'General_settings':{'threads':THREADS},
                         'Mappers':{'mapper_for_illumina':'boina',
                                    'mapper_for_sanger':'bwa'}}

        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                       configuration=configuration)
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        solexa = '@seq1\n'
        solexa += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq2\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq14\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq15\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq12\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq13\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq16\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq17\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTGCTAAGGGTTCTTGAGGATTTAT\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa += '@seq18\n'
        solexa += 'ATATGATTGAAGATATTTCTGGATGATCTAGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'

        sanger = '>seq3\n'
        sanger += 'GATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATAGGAGATACTGA'
        sanger += 'GATTCTGGAATCTCTGAGTTTCTGGGTTCAAGTTGCACTGACCATTGTTGGATTTGTAGA'
        sanger += 'TTGTTTCTTCTTTCATTAGGCATTGATTATGGGTAAATGCGTGGGTACATATAATATATA'
        sanger += 'TCTGTTGAATGCAATTTACACATTGACTGAGGAACAACATGAACATGGCAGCTTTCTCAA'
        sanger += 'AATTGAACCACAGAAGGCTTAAAAGCAAAGTCTTTGGAGAATCAGACTAAGCTTGAGA\n'
        sanger += '>seq4\n'
        sanger += 'TCCATACTTTACTCTATCTCTTTCTGTGTTTGGTAACACGCGAGGATTGGATGATAT'
        sanger += 'CTATATATTCTAATGTGGACTAAAAATGTGTGTGTGTGTATGAAGATGGGAAGCCGGAAG'
        sanger += 'TCATCAAGGAGAAAAGAGAGATAGACGACGAGAAGATGGCGTTGACGTTCAGAGGACTAG'
        sanger += 'AGGGTCATGTGATGGAGAAGTACAAGAAGTATGAGGTTATCTTACAGTTCATTCCCAAGT'
        sanger += 'CGAACGAAGGCTGCGTCTGCAAAGTCACTCTGATATGGGAGAATCGCAACGAAGACTCCC'
        sanger += '>seq5\n'
        sanger += 'TCCATACTTTACTCTATCTCTTTCTGTGTTTGGTAACACGCGAGGATTGGATGATAT'
        sanger += 'CTATATATTCTAAAGTGGACTAAAAATGTGTGTGTGTGTATGAAGATGGGAAGCCGGAAG'
        sanger += 'TCATCAAGGAGAAAAGAGAGATAGACGACGAGAAGATGGCGTTGACGTTCAGAGGACTAG'
        sanger += 'AGGGTCATGTGATGGAGAAGTACAAGAAGTATGAGGTTATCTTACAGTTCATTCCCAAGT'
        sanger += 'CGAACGAAGGCTGCGTCTGCAAAGTCACTCTGATATGGGAGAATCGCAACGAAGACTCCC'

        fpath_sanger = join(clean_reads_dir, 'lb_hola1.pl_sanger.sm_hola.fasta')
        fpath_solexa = join(clean_reads_dir,
                                    'lb_hola2.pl_illumina.sm_hola.sfastq')
        open(fpath_sanger, 'w').write(sanger)
        open(fpath_solexa, 'w').write(solexa)

        fpath_sanger2 = join(clean_reads_dir, 'lb_adios.pl_sanger.fasta')
        fpath_solexa2 = join(clean_reads_dir,
                                    'lb_adios.pl_illumina.sfastq')
        open(fpath_sanger2, 'w').write(sanger)
        open(fpath_solexa2, 'w').write(solexa)

        #the reference
        reference_dir = join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(join(DATA_DIR, 'blast/arabidopsis_genes')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping', silent=True)
        mapping_dir = join(project_dir, 'mapping')
        singular_mapping_dir = sorted(os.listdir(mapping_dir))[0]
        singular_mapping_dir = join(mapping_dir, singular_mapping_dir)
        assert exists(join(singular_mapping_dir, 'bams',
                            'by_readgroup', 'lb_hola2.pl_illumina.sm_hola.bam'))
        result_dir = join(mapping_dir, 'bams')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)
#        print result_dir_by_lib
#        raw_input()
        do_analysis(project_settings=settings_path, kind='merge_bams',
                    silent=True)
        assert exists(join(result_dir, 'merged.0.bam'))
        assert exists(join(result_dir, 'merged.0.bam.bai'))

        os.chdir('/tmp')
        test_dir.close()


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'TestBackboneMapping.test_mapping_analysis_boina']
    unittest.main()
