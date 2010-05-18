'''
This module is part of the ngs_backbone. It contains functions needed to run
all the analyses.

Created on 12/03/2010

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

import os, logging
from configobj import ConfigObj
from franklin.backbone.specifications import BACKBONE_DIRECTORIES

from franklin.backbone.annotation import DEFINITIONS as annot_defs
from franklin.backbone.cleaning import DEFINITIONS as clean_defs
from franklin.backbone.assembling import DEFINITIONS as assembly_defs
from franklin.backbone.mapping import DEFINITIONS as mapp_defs
from franklin.backbone.create_project import (validate_configuration,
                                              make_config_spec)

DEFINITIONS = [annot_defs, clean_defs, assembly_defs, mapp_defs]
BACKBONESPEC = ''

def get_analysis_especifications():
    'It groups all the especification'
    specifications = {}
    for spec in DEFINITIONS:
        for key, value in spec.items():
            specifications[key] = value
    return specifications

def _configure_logging(log_fpath, silent):
    'It prepares the logging infraestructure'
    logger = logging.getLogger('franklin')
    logger.setLevel(logging.INFO)
    #the format
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    #create a handler for the log file
    log_hand = logging.FileHandler(log_fpath)
    log_hand.setFormatter(formatter)
    logger.addHandler(log_hand)
    #create console handler and set level to info
    if not silent:
        console_hand = logging.StreamHandler()
        console_hand.setFormatter(formatter)
        logger.addHandler(console_hand)

def do_analysis(kind, project_settings=None, analysis_config=None, silent=False):
    'It does one of the predefined analyses'
    if project_settings is None:
        project_settings = os.path.join(os.getcwd(),
                                        BACKBONE_DIRECTORIES['config_file'])
        if not os.path.exists(project_settings):
            raise ValueError('Settings path not given and not found')

    if not analysis_config:
        analysis_config = {}

    config_spec = make_config_spec()
    settings = ConfigObj(project_settings, configspec=config_spec, unrepr=True)
    #settings = validate_configuration(settings)


    specifications = get_analysis_especifications()


    log_fpath = os.path.join(settings['General_settings']['project_path'],
                             BACKBONE_DIRECTORIES['log_file'])
    _configure_logging(log_fpath, silent)

    analysis_def = specifications[kind]

    analyzer_klass = analysis_def['analyzer']
    analyzer = analyzer_klass(project_settings=settings,
                        analysis_definition=analysis_def, silent=silent)

    analyzer.run()
