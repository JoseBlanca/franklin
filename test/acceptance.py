#!/usr/local/bin/python2.6
'It test the backbone whole process'

from optparse import OptionParser
import shutil, logging, os
from os.path import join, exists
from franklin.backbone.backbone_runner import do_analysis
from franklin.utils.misc_utils import NamedTemporaryDir, DATA_DIR

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
    repository_dir = join(DATA_DIR, 'acceptance')
    settings_path = prepare_conf(project_dir, repository_dir)
    choice = analysis
    if choice in ('cleaning', None):
        original_reads = join(project_dir,'reads/original')
        if exists(original_reads):
            os.remove(original_reads)
        reads = join(project_dir,'reads')
        if not exists(reads):
            os.mkdir(reads)
        os.symlink(join(repository_dir, 'cleaning'),
                   join(project_dir, 'reads/original'))
        analyses = ['clean_reads', 'clean_read_stats']
        run_analysis(analyses, settings_path)

    if choice in ('assembling', None):
        clean_reads_dir = join(project_dir, 'reads', 'cleaned')
        if os.path.exists(clean_reads_dir):
            shutil.rmtree(join(project_dir, 'reads'))
        os.mkdir(join(project_dir,'reads'))
        os.symlink(join(repository_dir, 'assembling'),
                   join(project_dir, 'reads/cleaned'))

        analyses = [ 'prepare_mira_assembly', 'mira_assembly',
                     'select_last_assembly']
        run_analysis(analyses, settings_path)

    if choice in ('mapping', None):
        clean_reads_dir = join(project_dir, 'reads', 'cleaned')
        if os.path.exists(clean_reads_dir):
            shutil.rmtree(join(project_dir, 'reads'))
        os.mkdir(join(project_dir, 'reads'))
        os.symlink(join(repository_dir, 'assembling'),
                   join(project_dir, 'reads/cleaned'))
        if exists(join(project_dir, 'mapping')):
            shutil.rmtree(join(project_dir, 'mapping'))
        os.makedirs(join(project_dir, 'mapping', 'reference'))
        os.symlink(join(repository_dir, 'mapping', 'reference.fasta'),
                   join(project_dir, 'mapping', 'reference', 'reference.fasta'))

        analyses = ['mapping', 'select_last_mapping', 'merge_bam',
                    'realign_bam']
        run_analysis(analyses, settings_path)

    if choice in ('snvs', None):
        annot_dir = join(project_dir, 'annotations')
        create_dir(annot_dir)
        annot_res = join(annot_dir, 'repr')
        os.mkdir(join(annot_dir, 'input'))
        os.mkdir(annot_res)
        os.symlink(join(repository_dir, 'snvs', 'reference.fasta'),
                   join(annot_dir, 'input', 'reference.fasta'))

        mapping_dir = join(project_dir, 'mapping')
        create_dir(mapping_dir)
        os.mkdir(join(mapping_dir, 'reference'))
        os.symlink(join(repository_dir, 'snvs', 'merged.bam'),
                   join(project_dir, 'mapping', 'merged.bam'))
        os.symlink(join(repository_dir, 'snvs', 'reference.fasta'),
                   join(project_dir, 'mapping', 'reference', 'reference.fasta'))
        analyses = ['annotate_snv', 'filter_snvs']
        run_analysis(analyses, settings_path)

    if choice in ('annotation', None):
        annot_dir = join(project_dir, 'annotations')
        if exists(join(annot_dir)):
            shutil.rmtree(annot_dir)
        os.mkdir(annot_dir)
        os.symlink(join(repository_dir, 'annotation', 'input'),
                   join(annot_dir, 'input'))
        os.symlink(join(repository_dir, 'annotation', 'blast'),
                   join(annot_dir, 'blast'))

        analyses = ['annotate_orf', 'annotate_microsatellite',
                    'annotate_go', 'annotate_description',
                    'annotate_orthologs', 'annotate_introns',
                    'write_annotation']
        run_analysis(analyses, settings_path)

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
    'It prepares the franklin configuration file'

    univec_database = join(DATA_DIR, 'blast', 'univec')
    estscan_matrix  = join(repository_dir, 'config_data', 'At.smat')
    tair7_seq       = join(repository_dir, 'annotation', 'tair7_genomic.fasta')

    out_fhand = open(join(project_dir, 'franklin.conf'), 'w')
    for line in open(join(repository_dir, 'franklin.conf')):
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
    os.symlink(join(repository_dir, 'config_data'),
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

