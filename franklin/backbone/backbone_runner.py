'''
Created on 12/03/2010

@author: peio
'''

import os, logging
from configobj import ConfigObj
from franklin.backbone.specifications import BACKBONE_DIRECTORIES

from franklin.backbone.annotation import DEFINITIONS as annot_defs
from franklin.backbone.cleaning import DEFINITIONS as clean_defs
from franklin.backbone.assembling import DEFINITIONS as assembly_defs
from franklin.backbone.mapping import DEFINITIONS as mapp_defs

DEFINITIONS = [annot_defs, clean_defs, assembly_defs, mapp_defs]

def get_analysis_especifications():
    'It groups all the especification'
    specifications = {}
    for spec in DEFINITIONS:
        for key, value in spec.items():
            specifications[key] = value
    return specifications

def _configure_logging(log_fpath):
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
    console_hand = logging.StreamHandler()
    console_hand.setFormatter(formatter)
    logger.addHandler(console_hand)

def do_analysis(kind, project_settings=None, analysis_config=None):
    'It does one of the predefined analyses'
    if project_settings is None:
        project_settings = os.path.join(os.getcwd(),
                                        BACKBONE_DIRECTORIES['config_file'])
        if not os.path.exists(project_settings):
            raise ValueError('Settings path not given and not found')

    if not analysis_config:
        analysis_config = {}

    settings = ConfigObj(project_settings, unrepr=True)
    specifications = get_analysis_especifications()
    log_fpath = os.path.join(settings['General_settings']['project_path'],
                             'franklin.log')
    _configure_logging(log_fpath)
    analysis_def = specifications[kind]

    analyzer_klass = analysis_def['analyzer']
    analyzer = analyzer_klass(project_settings=settings,
                        analysis_definition=analysis_def)

    analyzer.run()
