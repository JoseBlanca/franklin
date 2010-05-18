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
from configobj import ConfigObj, flatten_errors
from validate import (Validator, ValidateError, VdtTypeError, VdtValueError,
                      VdtValueTooSmallError, VdtValueTooBigError)
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

    config = make_config(project_path, name, config_data)
    config.filename = os.path.join(settings_path)

    #overwrite with the configuration given
    for section, config_info in configuration.items():
        for key, value in config_info.items():
            if not section in config:
                config[section] = {}
            config[section][key] = value

    # Validate Config
    validate_configuration(config)
    config.write()

    return settings_path

def validate_configuration(config, copy=True):
    'It validates the configuration'
    validator = Validator()
    #add function to validate some options
    validator.functions['integer_none_or_bool'] = is_integer_none_or_bool
    validator.functions['integer_or_none'] = is_integer_or_none
    validator.functions['string_or_none'] = is_string_or_none
    validator.functions['list_or_none'] = is_list_or_none
    result = config.validate(validator, copy=copy)
    errors = []
    if result != True:
        for (section_list, key, _) in flatten_errors(config, result):
            if key is not None:
                error_msg = 'The "%s" key in the section "%s" failed validation' % (key, ', '.join(section_list))
            else:
                error_msg = 'The following section was missing:%s ' % ', '.join(section_list)
            errors.append(error_msg)
        raise ValidateError('Check following error:\n %s' % '\n'.join(errors))
    else:
        config.write

def is_string_or_none(value):
    'This function validates the value and looks if it is an string or a None'
    kind = type(value)
    if  kind not in (type('str'), type(None)):
        raise VdtTypeError(value)
    return value

def is_list_or_none(value, min=None, max=None):
    'This function validates the value and looks if it is an list or a None'
    kind = type(value)
    if  kind not in (type([1]), type(None)):
        raise VdtTypeError(value)
    return value

def is_integer_none_or_bool(value, min=None, max=None):
    '''This function validates the value and looks if it is an integer a None or
    True'''

    kind = type(value)
    if kind not in (type(1), type(None), type(True)):
        raise VdtTypeError(value)

    if type(value) == type(1):
        if min is not None and value < int(min):
            raise VdtValueTooSmallError(value)
        if max is not None and value >int(max):
            raise VdtValueTooBigError(value)
    return value

def is_integer_or_none(value, min=None, max=None):
    '''This function validates the value and looks if it is an integer a None'''
    kind = type(value)
    if kind not in (type(1), type(None)):
        raise VdtTypeError(value)
    if isinstance(value, int):
        if min is not None and value < int(min):
            raise VdtValueTooSmallError(value)
        if max is not None and value >int(max):
            raise VdtValueTooBigError(value)
    return value


def make_config(project_path, name, config_data):
    'It creates the configobj'

    config_spec = make_config_spec()

    config = ConfigObj(configspec=config_spec, unrepr=True)

    config['General_settings'] = {}
    config['General_settings']['tmpdir'] = os.path.join(project_path, 'tmp')
    config['General_settings']['project_name'] = name
    config['General_settings']['project_path'] = project_path

    config['Cleaning'] = {}
    lucy_settings = os.path.join(config_data, 'lucy', 'lucy.conf')
    config['Cleaning']['lucy_settings'] = lucy_settings

    config['blast'] = {}
    config['blast']['nr'] = {}
    config['blast']['nr']['path'] = "Path_to_nr database"
    config['blast']['nr']['species'] = 'all'
    config['blast']['nr']['kind'] = 'prot'
    comments = []
    comments.append('Add as much blast databases as you need. Here a example')
    config['blast'].comments = {'nr':comments}


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

    return config

def make_config_spec():
    'It creates the backbone_config_spec '

    cfg ='''
    ['General_settings']
        tmpdir       = string
        project_name = string
        project_path = string
        threads      = integer_none_or_bool(min=1, default=None)

    ['Other_settings']
        default_sanger_quality = integer(default=20)
        java_memory            = integer(default=1024)
        picard_path            = string_or_none(default=None)
        gatk_path              = string_or_none(default=None)

    # Sample size to make stats
    ['Read_stats']
        sampling_size = integer(min=0, max=3000, default=1000)

    ['Cleaning']

        adaptors_file_454 = string_or_none(default=None)
        adaptors_file_sanger = string_or_none(default=None)
        adaptors_file_illumina = string_or_none(default=None)

        short_adaptors_sanger = list_or_none(default=None)
        short_adaptors_454 = list_or_none(default=None)
        short_adaptors_illumina = list_or_none(default=None)

        vector_database = string(default='UniVec')
        # Minimun seq length before removing the seq in the cleaning step
        [['min_seq_length']]
            454      = integer(min=0, max=400 ,default=100)
            sanger   = integer(min=0, max=400, default=100)
            illumina = integer(min=0, default=22)

        # Nucleotides to remove for each edge in each platform
        [['edge_removal']]
            454_left       = integer_or_none(default=None)
            454_right      = integer_or_none(default=None)
            sanger_left    = integer_or_none(default=None)
            sanger_right   = integer_or_none(default=None)
            illumina_left  = integer_or_none(default=None)
            illumina_right = integer_or_none(default=None)

    # Mira configuration.
    ['Mira']
        job_options = string_list(default=list('denovo', 'est'))
        general_settings = string_list(default=list('-AS:sd=1'))
        454_settings     = string_list(default=list("-LR:mxti=no", "-CO:rodirs=5", "-AL:mrs=80", "-OUT:sssip=0:stsip=0"))
        sanger_settings  = string_list(default=list('-AS:epoq=no', '-AS:bdq=30', '-CO:rodirs=5', '-AL:mrs=80', '-OUT:sssip=0:stsip=0'))


    ['Mappers']
        mapper_for_454 = string(default=bwa)
        mapper_for_illumina = string(default=bwa)
        mapper_for_solid = string(default=bwa)
        mapper_for_sanger = string(default=bwa)

    ['Sam_processing']
        add_default_qualities = boolean(default=True)

    ['Annotation']
        [['description_annotation']]
            # List of databases to use form description annotation.
            # This databases must be previously defined in blast section
            description_databases = string_list(default=list('somedb'))
        [['ortholog_annotation']]
            # List ofdatabses to look form orthoglogs
            # This databases must be previously defined in blast section
            ortholog_databases =  string_list(default=list('somedb'))
        [['Cdna_intron_annotation']]
            # Path to the genomic seqs
            genomic_seqs = string(default='path_to_seqs')
            genomic_db = string_or_none(default =None)

        [['orf_annotation']]
            estscan_matrix = string(default='path to estscan matrix')
        [['go_annotation']]
            blast_database = string(default='nr')
            java_memory = integer(default=2048)
            create_dat_file = boolean(default=False)

    ['blast']
        # Add as many blast databases as you need. You only need to add the
        # path and the species of the databases
        [[__many__]]
            path    = string(default="Path to database")
            species = string(default='species')


    ['Snvs']
        [['edge_removal']]
            454_left       = integer_or_none(default=None)
            454_right      = integer_or_none(default=None)
            sanger_left    = integer_or_none(default=None)
            sanger_right   = integer_or_none(default=None)
            illumina_left  = integer_or_none(default=None)
            illumina_right = integer_or_none(default=None)


'''

    cfg = cfg.splitlines()
    config_spec = ConfigObj(cfg, list_values=False,_inspec=True)
    return config_spec


