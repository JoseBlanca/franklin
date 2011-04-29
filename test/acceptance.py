#!/usr/bin/env python
'It test the backbone whole process'

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

from optparse import OptionParser
import shutil, logging, os
from os.path import join, exists
from franklin.backbone.backbone_runner import do_analysis
from franklin.backbone.specifications import BACKBONE_DIRECTORIES
from franklin.utils.misc_utils import NamedTemporaryDir, TEST_DATA_DIR

def test_backbone(analysis=None, analysis_dir=None):
    '''It tests the backbone infrastructure.

    If no analysis is given it will run all of them.
    If no analysis_dir is given a temporary one will be used.
    '''
    logger = logging.getLogger('franklin')
    if analysis_dir:
        analysis_fhand = None
        analysis_fpath = analysis_dir
    else:
        analysis_fhand = NamedTemporaryDir()
        analysis_fpath = analysis_fhand.name

    project_dir = analysis_fpath
    repository_dir = join(TEST_DATA_DIR, 'acceptance')
    settings_path = prepare_conf(project_dir, repository_dir)
    choice = analysis
    #choice = 'mapping'
    if choice in ('cleaning', None):
        original_reads = join(project_dir, 'reads/raw')
        if exists(original_reads):
            os.remove(original_reads)
        reads = join(project_dir, 'reads')
        if not exists(reads):
            os.mkdir(reads)
        shutil.copytree(join(repository_dir, 'cleaning'),
                        join(project_dir, 'reads/raw'))
        analyses = ['clean_reads', 'read_stats']
        run_analysis(analyses, settings_path)

    if choice in ('assembling', None):
        clean_reads_dir = join(project_dir, 'reads', 'cleaned')
        if os.path.exists(clean_reads_dir):
            shutil.rmtree(join(project_dir, 'reads'))
        os.mkdir(join(project_dir, 'reads'))
        shutil.copytree(join(repository_dir, 'assembling'),
                        join(project_dir, 'reads/cleaned'))

        analyses = [ 'prepare_mira_assembly', 'mira_assembly']
        run_analysis(analyses, settings_path)

    if choice in ('mapping', None):
        clean_reads_dir = join(project_dir, 'reads', 'cleaned')
        if os.path.exists(clean_reads_dir):
            shutil.rmtree(join(project_dir, 'reads'))
        os.mkdir(join(project_dir, 'reads'))
        shutil.copytree(join(repository_dir, 'assembling'),
                        join(project_dir, 'reads/cleaned'))
        if exists(join(project_dir, 'mapping')):
            shutil.rmtree(join(project_dir, 'mapping'))
        os.makedirs(join(project_dir, 'mapping', 'reference'))
        shutil.copy(join(repository_dir, 'mapping', 'reference.fasta'),
                    join(project_dir, 'mapping', 'reference',
                         'reference.fasta'))

        analyses = ['mapping', 'merge_bams', 'realign_bam']
        run_analysis(analyses, settings_path)

    if choice in ('snvs', None):
        annot_dir = join(project_dir, 'annotations')
        create_dir(annot_dir)
        annot_res = join(annot_dir, 'repr')
        os.mkdir(join(annot_dir, 'input'))
        os.mkdir(annot_res)
        shutil.copy(join(repository_dir, 'snvs', 'reference.fasta'),
                    join(annot_dir, 'input', 'reference.fasta'))

        mapping_dir = join(project_dir, 'mapping')
        create_dir(mapping_dir)
        os.mkdir(join(mapping_dir, 'reference'))
        shutil.copy(join(repository_dir, 'snvs', 'merged.bam'),
                   join(project_dir, 'mapping', 'merged.bam'))
        shutil.copy(join(repository_dir, 'snvs', 'reference.fasta'),
                   join(project_dir, 'mapping', 'reference', 'reference.fasta'))
        analyses = ['annotate_snvs', 'filter_snvs', 'annotation_stats',
                    'write_annotations']
        run_analysis(analyses, settings_path)

        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                            'reference.txt')
        result = open(stats_fpath).read()
        expected = '''Sequences with SNVs: 26
SNVs found: 72
SNV types:
\tinsertion: 2
\tdeletion: 7
\tcomplex: 1
\ttransition: 47
\ttransversion: 15
SNV locations:
\tunknown: 72'''
        #print  result
        assert expected in result

    if choice in ('annotation', None):
        annot_dir = join(project_dir, 'annotations')
        if exists(join(annot_dir)):
            shutil.rmtree(annot_dir)
        os.mkdir(annot_dir)
        shutil.copytree(join(repository_dir, 'annotation', 'input'),
                        join(annot_dir, 'input'))
        shutil.copytree(join(repository_dir, 'annotation', 'blast'),
                        join(annot_dir, 'blast'))

        analyses = ['annotate_orfs', 'annotate_microsatellites',
                    'annotate_gos', 'annotate_descriptions',
                    'annotate_orthologs', 'annotate_introns',
                    'annotate_prot_change',
                    'write_annotations', 'annotation_stats']
        run_analysis(analyses, settings_path)

        stats_fpath = join(project_dir, 'annotations', 'features', 'stats',
                            'tair7_cdna.st_nucl.txt')
        result = open(stats_fpath).read()
        expected = '''Number of sequences: 4
Sequences with description: 4
Sequences with ORF: 4
Number of ORFs: 4
Sequences with intron: 2
Number of introns: 3'''
        assert expected in result

    if not analysis_dir:
        analysis_fhand.close()

def create_dir(dir_, delete=True):
    'It creates the directory and if it exists it deletes previously'
    if delete and exists(join(dir_)):
        shutil.rmtree(dir_)
    os.mkdir(dir_)

def run_analysis(analyses, settings_path):
    'It runs the analyses and removes the files if fails'
    for analysis in analyses:
        print 'Running analysis: %s' % analysis
        do_analysis(project_settings=settings_path, kind=analysis,
                    silent=True)
    print "Test OK"

def prepare_conf(project_dir, repository_dir):
    'It prepares the backbone configuration file'

    univec_database = join(TEST_DATA_DIR, 'blast', 'univec+')
    estscan_matrix = join(repository_dir, 'config_data', 'At.smat')
    tair7_seq = join(repository_dir, 'annotation', 'tair7_genomic.fasta')

    out_fhand = open(join(project_dir, BACKBONE_DIRECTORIES['config_file']),
                     'w')
    for line in open(join(repository_dir, 'backbone.conf')):
        line = line.replace('/home/jope/test_backbone', project_dir)
        if 'UniVec' in line:
            line = line.replace('UniVec', univec_database)
        elif 'At.smat' in line:
            line = line.replace('/path/to/At.smat', estscan_matrix)
        elif 'tair7_genomic.fasta' in line:
            line = line.replace('/path/to/tair7_genomic.fasta', tair7_seq)
        out_fhand.write(line)
    out_fhand.close()
    config_data = join(project_dir, 'config_data')
    if exists(config_data):
        os.remove(config_data)
    shutil.copytree(join(repository_dir, 'config_data'),
                    join(project_dir, config_data))
    return out_fhand.name

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-d', '--work_dir', dest='work_dir',
                    help='directory to do the analysis')
    parser.add_option('-a', '--analysis', dest='analysis',
                      help='analysis to test')
    options = parser.parse_args()[0]
    return {'analysis':options.analysis, 'analysis_dir':options.work_dir}

if __name__ == '__main__':
    options = parse_options()
    test_backbone(**options)
