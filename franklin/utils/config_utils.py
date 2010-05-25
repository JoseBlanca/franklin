'''
ConfigObj utilities for franklin
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

import tempfile, shutil

def pretyfy_config(fpath):
    'It adds tabs to the config'
    new_config = tempfile.NamedTemporaryFile(delete=False)

    section_depth = 0
    for line in open(fpath):
        line = line.strip()
        if line.startswith('['):
            section_depth = line.count('[')

        if line.startswith('['):
            line = '    ' * (section_depth  - 1) + line
        else:
            line = '    ' * (section_depth - 0) + line

        line = line + '\n'
        if line.startswith('[') and line.count('[') == 1:
            new_config.write('\n')
        new_config.write(line)

    new_config.flush()
    shutil.move(new_config.name, fpath)
    new_config.close()

STRING = (basestring,)
INTEGER = (int,)
BOOLEAN = (bool,)
FLOAT = (float,)
NUMBER = (float, int)
INTEGER_OR_BOOL = (int, bool)

def _list(value, parameter, type_):
    'It checks that the values in the list are the given type_'
    if value is None:
        return
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

def validate_config(config, validation_spec):
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
            #is there a validation rule?
            if key not in validation_spec:
                continue
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

def add_default_values(config, defaults):
    'It adds the default values to the configuration'

    for key, value in defaults.items():
        if key == '__many__':
            #this is only used for the validation, not for the default
            continue

        #are we dealing with a dict or with just a value?
        if 'items' in dir(value):
            if key not in config:
                config[key] = {}
            add_default_values(config[key], value)
        else:
            #are there a default value?
            if key in config:
                continue
            if len(value) > 1:
                default_value = value[1]
                config[key] = default_value

