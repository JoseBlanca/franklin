#!/usr/bin/env python
'''
This script is used to create a backbone project
'''

from franklin.backbone.create_project import create_project
from optparse import OptionParser

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-p', '--project_name', dest='project_name',
                      help='Project name')
    parser.add_option('-d', '--directory', dest='directory',
                      help='working directory')
    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.project_name is None:
        parser.error('project_name is mandatory')
    else:
        project_name = options.project_name

    return project_name, options.directory

def main():
    'The main part'
    project_name, directory = set_parameters()
    create_project(project_name, directory=directory)

if __name__ == '__main__':
    main()
