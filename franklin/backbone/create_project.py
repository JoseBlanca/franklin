'''
Created on 29/01/2010

@author: jose
'''

import os
from configobj import ConfigObj
from franklin.backbone.analysis import BACKBONE_DIRECTORIES

def create_project(name, directory=None):
    'It creates the files that define a project'
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

    config = ConfigObj()
    config.filename = os.path.join(project_path,
                                   BACKBONE_DIRECTORIES['config_file'])

    config['General_settings'] = {}
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

    config['Snvs'] = {}
    config['Snvs']['snv_pipeline'] = ''
    config.write()

    return settings_path
