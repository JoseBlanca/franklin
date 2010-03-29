'''
Created on 29/03/2010

@author: jose
'''
import unittest, os
from os.path import join, exists

from configobj import ConfigObj

from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR
from franklin.backbone.create_project import create_project
from franklin.backbone.backbone_runner import do_analysis

class TestBackboneMapping(unittest.TestCase):
    'It tests the backbone'

    @staticmethod
    def test_mapping_analysis():
        'We can map the reads'
        test_dir = NamedTemporaryDir()
        project_name = 'backbone'
        settings_path = create_project(directory=test_dir.name,
                                       name=project_name,
                                      configuration={
                                                    'Snvs':{'min_quality':20},
                                'Sam_processing':{'add_default_qualities':True},
                                })
        project_dir = join(test_dir.name, project_name)
        #setup the original reads
        reads_dir = join(project_dir, 'reads')
        clean_reads_dir = join(reads_dir, 'cleaned')
        os.mkdir(reads_dir)
        os.mkdir(clean_reads_dir)

        solexa =  '@seq1\n'
        solexa += 'TCATTGAAAGTTGAAACTGATAGTAGCAGAGTTTTTTCCTCTGTTTGG\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIIIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq2\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq14\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq15\n'
        solexa += 'ATATGATTGAAGATATTTCTGGGCTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq12\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq13\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'
        solexa =  '@seq16\n'
        solexa += 'ATATGATTGAAGATATTTCTGGACTTTAAGGGTTCTTGAGGATTTATA\n'
        solexa += '+\n'
        solexa += 'IIIIIIHIIIIIIIIIIIIIIIZIIUJUAUGJUUJUDFAOUDJOFSUD\n'

        sanger  = '>seq3\n'
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

        #the reference
        reference_dir = join(project_dir, 'mapping/reference')
        os.makedirs(reference_dir)
        reference_fpath = join(reference_dir, 'reference.fasta')
        out = open(reference_fpath, 'w')
        for line in open(join(DATA_DIR, 'blast/arabidopsis_genes')):
            out.write(line)

        do_analysis(project_settings=settings_path, kind='mapping')
        mapping_dir = join(project_dir, 'mapping')
        singular_mapping_dir = sorted(os.listdir(mapping_dir))[0]
        singular_mapping_dir = join(mapping_dir, singular_mapping_dir)
        assert exists(join(singular_mapping_dir, 'result',
                                    'by_readgroup', 'lb_hola.pl_illumina.bam'))
        do_analysis(project_settings=settings_path, kind='select_last_mapping')
        result_dir = join(mapping_dir, 'result')
        assert exists(result_dir)
        result_dir_by_lib = join(result_dir, 'by_readgroup')
        assert exists(result_dir_by_lib)

        do_analysis(project_settings=settings_path, kind='merge_bam')
        assert exists(join(result_dir, 'merged.bam'))

        #we realign the mapping using GATK
        do_analysis(project_settings=settings_path, kind='realign_bam')

        annot_input_dir = join(project_dir, 'annotations', 'input')
        os.makedirs(annot_input_dir)
        os.symlink(reference_fpath, join(annot_input_dir, 'reference.fasta'))
        do_analysis(project_settings=settings_path, kind='annotate_snv')
        repr_fpath = join(project_dir, 'annotations', 'repr', 'reference.repr')
        assert "type='snv'" in  open(repr_fpath).read()

        do_analysis(project_settings=settings_path, kind='filter_snvs')
        do_analysis(project_settings=settings_path, kind='write_annotation')
        vcf_fpath = join(project_dir, 'annotations', 'result', 'reference.vcf')
        vcf = open(vcf_fpath).read()
        assert 'vks' in vcf
        assert 'AT5G19860.1' in vcf

        os.chdir('/tmp')
        test_dir.close()

if __name__ == "__main__":
    #import sys;sys.argv = ['TestBackbone.test_mapping_analysis']#, 'Test.testName']
    unittest.main()
