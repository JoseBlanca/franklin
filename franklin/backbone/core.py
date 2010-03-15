'''
Created on 12/03/2010

@author: peio
'''

import os
from configobj import ConfigObj
from franklin.backbone.analysis import ANALYSIS_DEFINITIONS
from franklin.backbone.specifications import BACKBONE_DIRECTORIES
from franklin.backbone.annotation import SPECIFICATIONS


def get_analysis_especifications():
    'It groups all the especification'
    specifications = {}
    for spec in [SPECIFICATIONS, ANALYSIS_DEFINITIONS]:
        for key, value in spec.items():
            specifications[key] = value
    return specifications

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
    analysis_def = specifications[kind]

    analyzer_klass = analysis_def['analyzer']
    analyzer = analyzer_klass(project_settings=settings,
                        analysis_definition=analysis_def)

    analyzer.run()
