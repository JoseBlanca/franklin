'''
Created on 29/01/2010

@author: jose
'''

import os
from configobj import ConfigObj
from franklin.backbone.analysis import BACKBONE_DIRECTORIES

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
    settings_path = os.path.join(project_path, 'backbone.conf')
    config_data = os.path.join(project_path, 'config_data')

    config = ConfigObj(unrepr=True)
    config.filename = os.path.join(project_path,
                                   BACKBONE_DIRECTORIES['config_file'])

    config['General_settings'] = {}
    config['General_settings']['tmpdir'] = os.path.join(project_path, 'tmp')
    config['General_settings']['project_name'] = name
    config['General_settings']['project_path'] = project_path

    config['Cleaning'] = {}
    config['Cleaning']['vector_database'] = 'UniVec'
    comments = []
    comments.append('adaptors_file_454 = /some/adaptors/fasta/file')
    comments.append('adaptors_file_sanger = /some/adaptors/fasta/file')
    comments.append('adaptors_file_illumina = /some/adaptors/fasta/file')
    config['Cleaning'].comments = {'vector_database':comments}

    lucy_settings = os.path.join(config_data, 'lucy', 'lucy.conf')
    config['Cleaning']['lucy_settings'] = lucy_settings

    config['Mira'] = {}
    config['Mira']['job_options']     = ['denovo', 'est']
    config['Mira']['454_settings']    = ['-LR:mxti=no']
    config['Mira']['sanger_settings'] = ['-AS:epoq=no', '-AS:bdq=30']

    config['Mappers'] = {}
    config['Mappers']['mapper_for_454'] = 'bwa'
    config['Mappers']['mapper_for_illumina'] = 'bwa'
    config['Mappers']['mapper_for_solid'] = 'bwa'
    config['Mappers']['mapper_for_sanger'] = 'bwa'

    config['snv_filters'] = {}
    config['snv_filters']['filter1'] = {}
    config['snv_filters']['filter1']['name'] = 'uniq_contiguous'
    config['snv_filters']['filter1']['use']  = False

    config['snv_filters']['filter2'] = {}
    config['snv_filters']['filter2']['name']     = 'close_to_intron'
    config['snv_filters']['filter2']['use']      = False
    config['snv_filters']['filter2']['distance'] = 30

    config['snv_filters']['filter3'] = {}
    config['snv_filters']['filter3']['name']            = 'high_variable_region'
    config['snv_filters']['filter3']['use']             = False
    config['snv_filters']['filter3']['max_variability'] = 0.6
    config['snv_filters']['filter3']['window']          = None

    config['snv_filters']['filter4'] = {}
    config['snv_filters']['filter4']['name']      = 'close_to_snv'
    config['snv_filters']['filter4']['use']       = False
    config['snv_filters']['filter4']['proximity'] = 60

    config['snv_filters']['filter5'] = {}
    config['snv_filters']['filter5']['name']     = 'close_to_limit'
    config['snv_filters']['filter5']['use']      = False
    config['snv_filters']['filter5']['distance'] = 60

    config['snv_filters']['filter6'] = {}
    config['snv_filters']['filter6']['name']      = 'maf'
    config['snv_filters']['filter6']['use']       = False
    config['snv_filters']['filter6']['frecuency'] = 0.8

    config['snv_filters']['filter7'] = {}
    config['snv_filters']['filter7']['name'] = 'by_kind'
    config['snv_filters']['filter7']['use']  = True
    config['snv_filters']['filter7']['kind'] = 0 # snp

    config['snv_filters']['filter8'] = {}
    config['snv_filters']['filter8']['name']        = 'cap_enzyme'
    config['snv_filters']['filter8']['use']         = False
    config['snv_filters']['filter8']['all_enzymes'] = True

    config['snv_filters']['filter9'] = {}
    config['snv_filters']['filter9']['name']       = 'is_variable_in_rg'
    config['snv_filters']['filter9']['use']        = False
    config['snv_filters']['filter9']['group_kind'] = 'read_groups'
    config['snv_filters']['filter9']['groups']     = None

    config['snv_filters']['filter10'] = {}
    config['snv_filters']['filter10']['name']       = 'is_variable_in_lb'
    config['snv_filters']['filter10']['use']        = False
    config['snv_filters']['filter10']['group_kind'] = 'libraries'
    config['snv_filters']['filter10']['groups']     = None

    config['snv_filters']['filter11'] = {}
    config['snv_filters']['filter11']['name']       = 'is_variable_in_sm'
    config['snv_filters']['filter11']['use']        = False
    config['snv_filters']['filter11']['group_kind'] = 'samples'
    config['snv_filters']['filter11']['groups']     = None

    #overwrite with the configuration given
    outputs = ['vcf']
    for section, config_info in configuration.items():
        for key, value in config_info.items():
            if not section in config:
                config[section] = {}
            config[section][key] = value

    config.write()

    return settings_path
