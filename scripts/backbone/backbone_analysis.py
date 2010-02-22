#!/usr/bin/env python
'''
This script is used to run any of the backbone analysis.
'''

from franklin.backbone.analysis import do_analysis
from optparse import OptionParser
import os

ACTION_MSG = ''


def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-a', '--action', dest='action',  help=ACTION_MSG)
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

    if options.settings_fpath is None:
        if os.path.exists('settings.py'):
            settings_fpath = os.path.abspath('settings.py')
        else:
            parser.error('Settings file path is mandatory')
    else:
        settings_fpath = options.settings_fpath

    return action, settings_fpath

def main():
    'The main part'
    action, settings_fpath = set_parameters()
    do_analysis(project_settings=settings_fpath, kind=action)

if __name__ == '__main__':
    main()
