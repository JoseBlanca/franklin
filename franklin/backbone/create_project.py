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
from validate import Validator

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

    config = _make_default_config(project_path, name, config_data)
    config.filename = os.path.join(settings_path)

    #overwrite with the configuration given
    for section, config_info in configuration.items():
        for key, value in config_info.items():
            if not section in config:
                config[section] = {}
            config[section][key] = value

    config.write()

    print open(config.filename).read()

    return settings_path

def create_configuration(config_path):
    'It returns a validated configObj'

    config = ConfigObj(config_path, unrepr=True)

    _add_default_values(config, DEFAULT_CONFIGURATION)

    _validate_config(config, DEFAULT_CONFIGURATION)

    return config

def _make_default_config(project_path, name, config_data):
    'It creates the configobj adding some default values'

    config = ConfigObj(unrepr=True)

    _add_default_values(config, DEFAULT_CONFIGURATION)

    _add_some_settings(config, project_path, name, config_data)

    _validate_config(config, DEFAULT_CONFIGURATION)

    return config

STRING = (basestring,)
INTEGER = (int,)
BOOLEAN = (bool,)
FLOAT = (float,)
NUMBER = (float, int)
INTEGER_OR_BOOL = (int, bool)

def _list(value, parameter, type_):
    'It checks that the values in the list are the given type_'
    if not isinstance(value, list) and not isinstance(value, tuple):
        msg = 'In parameter %s, value %s is not a list' % (parameter,
                                                            str(value))
        raise ValueError(msg)
    for item in value:
        _check_type(item, parameter, type_)

STRING_LIST = (_list, STRING)
FLOAT_LIST = (_list, FLOAT)
NUMBER_LIST = (_list, NUMBER)

def _check_type(value, parameter, types):
    'It check that the value is an string or None'
    if value is None:
        return
    if not any([isinstance(value, type_) for type_ in types]):
        msg = 'In parameter %s, value %s is not a %s' % (parameter,
                                                         str(value),
                                                         str(type_))
        raise ValueError(msg)

def _validate_config(config, validation_spec):
    'It raises an error if some value is not ok'
    errors = _validate_config_recursive(config, validation_spec)
    if errors:
        msg = 'Erros found in the configuration for the following '
        msg += 'sections:values:type -> '
        errors = ', '.join(['%s:%s:%s' % (error) for error in errors])
        msg += errors
        raise ValueError(msg)

def _validate_config_recursive(config, validation_spec, errors=None):
    'It returns the errors found in the configuration'
    validator_strs = {(basestring,): 'string',
                      (int,): 'integer',
                      (bool,): 'bolean',
                      (float,): 'float',
                      (float, int): 'number',
                      (int, bool): 'integer or boolean',
                      (_list, STRING): 'list of strings',
                      (_list, FLOAT): 'list of floats',
                      (_list, NUMBER): 'list of numbers'}

    if errors is None:
        errors = []
    for key, value in config.items():
        #if key == '__many__':
            #this is only used for the validation, not for the default
            #continue

        #are we dealing with a dict or with just a value?
        if 'items' in dir(value):
            if '__many__' in validation_spec:
                spec_value = validation_spec['__many__']
            elif key not in validation_spec:
                #some sections do not require validation
                spec_value = None
            else:
                spec_value = validation_spec[key]
            if spec_value is not None:
                _validate_config_recursive(value, spec_value, errors)
        else:
            #are there a default value?
            if key == 'kind':
                pass
            spec_value = validation_spec[key]
            validator = validation_spec[key][0]
            validator_str = validator_strs[validator]
            try:
                if validator in (STRING, INTEGER, BOOLEAN, FLOAT,
                                 INTEGER_OR_BOOL):
                    _check_type(value, key, validator)
                elif validator in (STRING_LIST, FLOAT_LIST, NUMBER_LIST):
                    validator[0](value, key, validator[1])
                else:
                    raise RuntimeError('Unknown validator ' + str(validator))

            except ValueError:
                errors.append((key, value, validator_str))
    return errors

def _add_default_values(config, defaults):
    'It adds the default values to the configuration'

    for key, value in defaults.items():
        if key == '__many__':
            #this is only used for the validation, not for the default
            continue

        #are we dealing with a dict or with just a value?
        if 'items' in dir(value):
            if key not in config:
                config[key] = {}
            _add_default_values(config[key], value)
        else:
            #are there a default value?
            if key in config:
                continue
            if len(value) > 1:
                default_value = value[1]
                config[key] = default_value

def _add_some_settings(config, project_path, name, config_data):
    'It adds some extra default parameters'

    config['General_settings'] = config['General_settings']
    config['General_settings']['project_name'] = name
    config['General_settings']['project_path'] = project_path
    config['General_settings']['tmpdir'] = os.path.join(project_path, 'tmp')

    config['blast'] = config['blast']
    config['blast']['nr'] = {}
    config['blast']['nr']['path'] = "/path/to/nr/database"
    config['blast']['nr']['species'] = 'all'
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

from franklin.utils.misc_utils import OrderedDict

DEFAULT_CONFIGURATION = OrderedDict([
           ('General_settings',
                OrderedDict([
                    ('tmpdir', (STRING, None)),
                    ('project_name', (STRING, None)),
                    ('project_path', (STRING, None)),
                    ('threads', (INTEGER_OR_BOOL, None)),
                ]),
            ),
           ('Other_settings',
                OrderedDict([
                    ('default_sanger_quality', (INTEGER, 20)),
                    ('java_memory', (INTEGER, None)),
                    ('picard_path', (STRING, None)),
                    ('gatk_path',   (STRING, None)),
                ])
            ),
           ('Cleaning',
                OrderedDict([
                    ('adaptors_file_454', (STRING, None)),
                    ('adaptors_file_sanger', (STRING, None)),
                    ('adaptors_file_illumina', (STRING, None)),
                    ('short_adaptors_sanger', (STRING_LIST, [])),
                    ('short_adaptors_454', (STRING_LIST, [])),
                    ('short_adaptors_illumina', (STRING_LIST, [])),
                    ('vector_database', (STRING, 'UniVec')),
                    ('min_seq_length',{
                                      '454' : (INTEGER, 100),
                                      'sanger': (INTEGER, 100),
                                      'illumina': (INTEGER, 22),
                    }),
                    ('edge_removal',
                         OrderedDict([
                            ('454_left', (INTEGER, None)),
                            ('454_right', (INTEGER, None)),
                            ('sanger_left', (INTEGER, None)),
                            ('sanger_right', (INTEGER, None)),
                            ('illumina_left', (INTEGER, None)),
                            ('illumina_right', (INTEGER, None)),
                        ]),
                    ),
                    ('lucy', {
                             'vector_settings': (STRING, None),
                             'bracket': (NUMBER_LIST, [10, 0.02]),
                             'window': (NUMBER_LIST, [50, 0.08, 10, 0.3]),
                             'error': (NUMBER_LIST, [0.015, 0.015])
                    })
               ]),
            ),
           ('Mira',
                OrderedDict([
                   ('job_options', (STRING_LIST, ['denovo', 'est'])),
                   ('general_settings', (STRING_LIST, ['-AS:sd=1'])),
                   ('454_settings', (STRING_LIST, ["-LR:mxti=no",
                                                   "-CO:rodirs=5",
                                                   "-AL:mrs=80",
                                                   "-OUT:sssip=0:stsip=0"])),
                   ('sanger_settings', (STRING_LIST, ['-AS:epoq=no',
                                                      '-AS:bdq=30',
                                                      '-CO:rodirs=5',
                                                      '-AL:mrs=80',
                                                      '-OUT:sssip=0:stsip=0'])),
                ]),
           ),
           ('Mappers', {
                'mapper_for_454': (STRING, 'bwa'),
                'mapper_for_illumina': (STRING, 'bwa'),
                'mapper_for_solid': (STRING, 'bwa'),
                'mapper_for_sanger': (STRING, 'bwa'),
           }),
           ('Sam_processing',{
                'add_default_qualities': (BOOLEAN, True),
            }),
           ('Sam_stats',{
                'sampling_size': (INTEGER, None),
            }),
          ('Annotation', {
                'description_annotation':{
                    'description_databases': (STRING_LIST, [])},
                'ortholog_annotation': {
                    'ortholog_databases': (STRING_LIST, [])},
                'Cdna_intron_annotation': {
                    'genomic_seqs': (STRING, None),
                    'genomic_db' : (STRING, None),},
                'orf_annotation':{
                    'estscan_matrix' : (STRING, 'path to estscan matrix')},
                'go_annotation':{
                    'blast_database': (STRING, 'nr'),
                    'java_memory': (INTEGER, 2048),
                    'create_dat_file': (BOOLEAN, False),
                    'prop_fpath': (STRING, None)}
                },
            ),
        ('blast', {
            '__many__':{
                'path': (STRING, "/Path/to/database"),
                'species': (STRING, 'species_name')},
            }),
        ('Snvs',
            OrderedDict([
                ('min_quality', (INTEGER, 45)),
                ('min_mapq', (INTEGER, 15)),
                ('min_num_alleles', (INTEGER, 1)),
                ('edge_removal',
                         OrderedDict([
                            ('454_left', (INTEGER, None)),
                            ('454_right', (INTEGER, None)),
                            ('sanger_left', (INTEGER, None)),
                            ('sanger_right', (INTEGER, None)),
                            ('illumina_left', (INTEGER, None)),
                            ('illumina_right', (INTEGER, None)),
                        ])
                )
            ])
        )
    ])
