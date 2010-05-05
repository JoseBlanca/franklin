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

import unittest, os
from os.path import join, exists

from configobj import ConfigObj

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis

THREADS = 2

class TestBackboneMapping(unittest.TestCase):
    'It tests the backbone'

    @staticmethod
    def test_mapping_analysis():
        'We can map the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        blastdb_seq = os.path.join(DATA_DIR, 'blast', 'arabidopsis_genes')
        snv_filters = {'filter1':{'name':'uniq_contiguous', 'use':True,
                                  'genomic_db':blastdb_seq,
                                  'genomic_seqs_fpath':blastdb_seq},
                       'filter12':{'name':'ref_not_in_list', 'use':True,
                                'list_path':os.path.join(DATA_DIR, 'cos_list')},
                       'filter10':{'name': 'variable_in_sm',
                                   'step_name': 'is_variable', 'use':True,
                                   'group_kind':'libraries',
                                   'groups':['hola']}}

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

        fpath_sanger = join(clean_reads_dir, 'lb_hola.pl_sanger.fasta')
        fpath_solexa = join(clean_reads_dir,
                                    'lb_hola.pl_illumina.sfastq')
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
        assert exists(join(singular_mapping_dir, 'result',
                                    'by_readgroup', 'lb_hola.pl_illumina.bam'))
        do_analysis(project_settings=settings_path, kind='select_last_mapping',
                    silent=True)
        result_dir = join(mapping_dir, 'result')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)
        f1 = open(join(result_dir, 'merged.bam'), 'w')
        f1.close()
        do_analysis(project_settings=settings_path, kind='merge_bam',
                    silent=True)
        assert exists(join(result_dir, 'merged.0.bam'))


        #we realign the mapping using GATK
        do_analysis(project_settings=settings_path, kind='realign_bam',
                    silent=True)

        assert exists(join(result_dir, 'merged.1.bam'))
        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)
        os.symlink(reference_fpath, join(annot_input_dir, 'reference.fasta'))
        do_analysis(project_settings=settings_path, kind='annotate_snv',
                    silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr',
                          'reference.0.repr')
        assert "type='snv'" in  open(repr_fpath).read()

        do_analysis(project_settings=settings_path, kind='filter_snvs',
                    silent=True)
        repr_fpath = join(project_dir, 'annotations', 'repr',
                          'reference.1.repr')
        result = open(repr_fpath).read()
        assert "type='snv'" in result

        do_analysis(project_settings=settings_path, kind='write_annotation',
                    silent=True)
        vcf_fpath = join(project_dir, 'annotations', 'result',
                         'reference.vcf')
        vcf = open(vcf_fpath).read()
        assert 'VKS' in vcf
        assert 'AT5G19860.1' in vcf

        os.chdir('/tmp')
        test_dir.close()

if __name__ == "__main__":
    #import sys;sys.argv = ['TestBackbone.test_mapping_analysis']#, 'Test.testName']
    unittest.main()
