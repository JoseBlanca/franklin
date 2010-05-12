'''
This module is part of ngs_backbone. It provides of function to create new
backbone projects.

Created on 29/01/2010

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

import os
from configobj import ConfigObj
from franklin.backbone.specifications import BACKBONE_DIRECTORIES

def create_project(name, directory=None, configuration=None):
    'It creates the files that define a project'
    if configuration is None:
        configuration = {}
    if not directory:
        directory = os.getcwd()
    project_path = os.path.join(os.path.abspath(directory), name)
    if os.path.exists(project_path):
        raise ValueError('Directory already exists: ' + project_path)

    #create the directory
    os.mkdir(project_path)

    #create the settings
    settings_path = os.path.join(project_path,
                                 BACKBONE_DIRECTORIES['config_file'])
    config_data = os.path.join(project_path,
                               BACKBONE_DIRECTORIES['external_software_config'])

    config = ConfigObj(unrepr=True)
    config.filename = os.path.join(settings_path)

    config['General_settings'] = {}
    config['General_settings']['tmpdir'] = os.path.join(project_path, 'tmp')
    config['General_settings']['project_name'] = name
    config['General_settings']['project_path'] = project_path
    config['General_settings']['num_threads'] = None

    config['Other_settings'] = {}
    config['Other_settings']['default_sanger_quality'] = 20
    config['Other_settings']['java_memory'] = 2048

    config['Read_stats'] = {}
    config['Read_stats']['sampling_size'] = 1000

    config['Cleaning'] = {}

    config['Cleaning']['vector_database'] = 'UniVec'
    comments = []
    comments.append('adaptors_file_454 = /some/adaptors/fasta/file')
    comments.append('adaptors_file_sanger = /some/adaptors/fasta/file')
    comments.append('adaptors_file_illumina = /some/adaptors/fasta/file')

    comments.append('words_to_remove_sanger = [someword, another_word]')
    comments.append('words_to_remove_454 = [someword, another_word]')
    comments.append('words_to_remove_illumina = [someword, another_word]')
    config['Cleaning'].comments = {'vector_database':comments}

    config['Cleaning']['min_seq_length'] = {}
    config['Cleaning']['min_seq_length']['454'] = 100
    config['Cleaning']['min_seq_length']['sanger'] = 100
    config['Cleaning']['min_seq_length']['illumina'] = 22

    config['Cleaning']['edge_removal'] = {}
    config['Cleaning']['edge_removal']['454_left'] = None
    config['Cleaning']['edge_removal']['454_right'] = None
    config['Cleaning']['edge_removal']['sanger_left'] = None
    config['Cleaning']['edge_removal']['sanger_right'] = None
    config['Cleaning']['edge_removal']['illumina_left'] = None
    config['Cleaning']['edge_removal']['illumina_right'] = None

    lucy_settings = os.path.join(config_data, 'lucy', 'lucy.conf')
    config['Cleaning']['lucy_settings'] = lucy_settings

    config['Mira'] = {}
    config['Mira']['job_options'] = ['denovo', 'est']
    config['Mira']['general_settings'] = ['-AS:sd=1']
    config['Mira']['454_settings'] = ['-LR:mxti=no', '-CO:rodirs=5',
                                         '-AL:mrs=80',
                                         '-OUT:sssip=0:stsip=0']
    config['Mira']['sanger_settings'] = ['-AS:epoq=no', '-AS:bdq=30',
                                         '-CO:rodirs=5', '-AL:mrs=80',
                                         '-OUT:sssip=0:stsip=0']

    config['Mappers'] = {}
    config['Mappers']['mapper_for_454'] = 'bwa'
    config['Mappers']['mapper_for_illumina'] = 'bwa'
    config['Mappers']['mapper_for_solid'] = 'bwa'
    config['Mappers']['mapper_for_sanger'] = 'bwa'

    config['Sam_processing'] = {}
    config['Sam_processing']['add_default_qualities'] = True

    config['Annotation'] = {}
    config['Annotation']['description_annotation'] = {}
    config['Annotation']['description_annotation']['description_databases'] = ['somedb']
    config['Annotation']['ortholog_annotation'] = {}
    config['Annotation']['ortholog_annotation']['ortholog_databases'] = ['somedb']
    config['Annotation']['Cdna_intron_annotation'] = {}
    config['Annotation']['Cdna_intron_annotation']['genomic_db'] = 'path_to_blastdb'
    config['Annotation']['Cdna_intron_annotation']['genomic_seqs'] = 'path_to_seqs'
    config['Annotation']['orf_annotation'] = {}
    config['Annotation']['orf_annotation']['estscan_matrix'] = 'path to estscan matrix'
    config['Annotation']['go_annotation'] = {}
    config['Annotation']['go_annotation']['blast_database'] = 'nr'
    config['Annotation']['go_annotation']['java_memory'] = 2048
    config['Annotation']['go_annotation']['create_dat_file'] = False


    config['blast'] = {}
    config['blast']['nr'] = {}
    config['blast']['nr']['path'] = "Path_to_nr database"
    config['blast']['nr']['species'] = 'all'
    config['blast']['nr']['kind'] = 'prot'
    comments = []
    comments.append('Add as much blast databases as you need. Here a example')
    config['blast'].comments = {'nr':comments}

    config['Snvs'] = {}
    config['Snvs']['edge_removal'] = {}
    config['Snvs']['edge_removal']['454_left'] = None
    config['Snvs']['edge_removal']['454_right'] = None
    config['Snvs']['edge_removal']['sanger_left'] = None
    config['Snvs']['edge_removal']['sanger_right'] = None
    config['Snvs']['edge_removal']['illumina_left'] = None
    config['Snvs']['edge_removal']['illumina_right'] = None

    config['snv_filters'] = {}
    config['snv_filters']['filter1'] = {}
    config['snv_filters']['filter1']['name'] = 'uniq_contiguous'
    config['snv_filters']['filter1']['use'] = False
    config['snv_filters']['filter1']['genomic_db'] = 'path to blast db'
    config['snv_filters']['filter1']['genomic_seqs_fpath'] = 'path to seqs file'

    config['snv_filters']['filter2'] = {}
    config['snv_filters']['filter2']['name'] = 'close_to_intron'
    config['snv_filters']['filter2']['use'] = False
    config['snv_filters']['filter2']['distance'] = 30

    config['snv_filters']['filter3'] = {}
    config['snv_filters']['filter3']['name'] = 'high_variable_region'
    config['snv_filters']['filter3']['use'] = False
    config['snv_filters']['filter3']['max_variability'] = 0.06
    config['snv_filters']['filter3']['window'] = None

    config['snv_filters']['filter4'] = {}
    config['snv_filters']['filter4']['name'] = 'close_to_snv'
    config['snv_filters']['filter4']['use'] = False
    config['snv_filters']['filter4']['distance'] = 60

    config['snv_filters']['filter5'] = {}
    config['snv_filters']['filter5']['name'] = 'close_to_limit'
    config['snv_filters']['filter5']['use'] = False
    config['snv_filters']['filter5']['distance'] = 60

    config['snv_filters']['filter6'] = {}
    config['snv_filters']['filter6']['name'] = 'maf'
    config['snv_filters']['filter6']['use'] = False
    config['snv_filters']['filter6']['frequency'] = 0.8
    config['snv_filters']['filter6']['group_kind'] = 'read_groups'
    config['snv_filters']['filter6']['groups'] = []

    config['snv_filters']['filter7'] = {}
    config['snv_filters']['filter7']['name'] = 'by_kind'
    config['snv_filters']['filter7']['use'] = True
    config['snv_filters']['filter7']['kind'] = 0 # snp

    config['snv_filters']['filter8'] = {}
    config['snv_filters']['filter8']['name'] = 'cap_enzyme'
    config['snv_filters']['filter8']['use'] = False
    config['snv_filters']['filter8']['all_enzymes'] = True

    config['snv_filters']['filter9'] = {}
    config['snv_filters']['filter9']['name'] = 'is_variable_in_rg'
    config['snv_filters']['filter9']['step_name'] = 'is_variable'
    config['snv_filters']['filter9']['use'] = False
    config['snv_filters']['filter9']['group_kind'] = 'read_groups'
    config['snv_filters']['filter9']['groups'] = []

    config['snv_filters']['filter10'] = {}
    config['snv_filters']['filter10']['name'] = 'is_variable_in_lb'
    config['snv_filters']['filter10']['step_name'] = 'is_variable'
    config['snv_filters']['filter10']['use'] = False
    config['snv_filters']['filter10']['group_kind'] = 'libraries'
    config['snv_filters']['filter10']['groups'] = []

    config['snv_filters']['filter11'] = {}
    config['snv_filters']['filter11']['name'] = 'is_variable_in_sm'
    config['snv_filters']['filter11']['step_name'] = 'is_variable'
    config['snv_filters']['filter11']['use'] = False
    config['snv_filters']['filter11']['group_kind'] = 'samples'
    config['snv_filters']['filter11']['groups'] = []

    config['snv_filters']['filter12'] = {}
    config['snv_filters']['filter12']['name'] = 'ref_not_in_list'
    config['snv_filters']['filter12']['use'] = False
    config['snv_filters']['filter12']['list_path'] = 'path_to_file_with_list'

    #overwrite with the configuration given
    outputs = ['vcf']
    for section, config_info in configuration.items():
        for key, value in config_info.items():
            if not section in config:
                config[section] = {}
            config[section][key] = value

    config.write()

    return settings_path
