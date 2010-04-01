#!/usr/bin/env python
'''
This script is used to run any of the backbone analysis.
'''

import os, logging, cgitb, datetime

from franklin.backbone.backbone_runner import (do_analysis,
                                               get_analysis_especifications)
from franklin.backbone.specifications import BACKBONE_DIRECTORIES
from optparse import OptionParser

ERROR_DIR = os.path.join(BACKBONE_DIRECTORIES['error_dir'])
cgitb.enable(display=0, format='text', logdir=ERROR_DIR)

def _get_available_analyses():
    'It return available analyses'
    analyses = get_analysis_especifications().keys()
    return ", ".join(analyses)

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    msg = 'Posible analyses to run: %s' % _get_available_analyses()
    parser.add_option('-a', '--action', dest='action',  help=msg)
    parser.add_option('-s', '--settings', dest='settings',
                      help='Settings file path')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.action is None:
        parser.error('Action is mandatory')
    else:
        action = options.action

    if options.settings is None:
        config_fname = BACKBONE_DIRECTORIES['config_file']
        if os.path.exists(config_fname):
            settings_fpath = os.path.abspath(config_fname)
        else:
            msg = 'Settings file path is mandatory if %s is not found' % \
                                                                    config_fname
            parser.error(msg)
    else:
        settings_fpath = options.settings

    return action, settings_fpath

def main():
    'The main part'
    action, settings_fpath = set_parameters()
    logger = logging.getLogger('franklin')
    try:
        start_time = datetime.datetime.today()
        do_analysis(project_settings=settings_fpath, kind=action)
        time_elapsed = datetime.datetime.today() - start_time
        logger.info('Time elapsed %s' % str(time_elapsed))
    except Exception as error:
        logger.exception(error)
        if not os.path.exists(ERROR_DIR):
            os.mkdir(ERROR_DIR)
        cgitb.handler()
        raise

if __name__ == '__main__':
    main()
