'''
It test backbone_runner module

Created on 12/05/2010

@author: peio
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

import unittest
from configobj import ConfigObj, flatten_errors
from tempfile import NamedTemporaryFile
from validate import Validator
from franklin.backbone.create_project import make_config, validate_configuration

class TestValidator(unittest.TestCase):
    'It test configuration validator'

    @staticmethod
    def test_backbone_validator():
        'It test that we validate the config'
        project_path, name, config_data = 'path', 'test', 'path/config'
        config = make_config(project_path, name, config_data)
        validate_configuration(config, copy=True)
        assert config['General_settings']['project_name'] == 'test'
        assert config['Other_settings']['java_memory'] == 1024
        fhand = NamedTemporaryFile()
        fpath = fhand.name
        config.filename = fpath
        config.write()
        config_content = open(fhand.name).read()
        assert   'Nucleotides' in  config_content


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
