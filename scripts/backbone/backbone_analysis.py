#!/usr/bin/env python
'''
This script is used to run any of the backbone analysis.
'''

from franklin.backbone.backbone_runner import (do_analysis,
                                               get_analysis_especifications)

from optparse import OptionParser
import os


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
        if os.path.exists('backbone.conf'):
            settings_fpath = os.path.abspath('backbone.conf')
        else:
            parser.error('Settings file path is mandatory')
    else:
        settings_fpath = options.settings

    return action, settings_fpath

def main():
    'The main part'
    action, settings_fpath = set_parameters()
    do_analysis(project_settings=settings_fpath, kind=action)

if __name__ == '__main__':
    main()
