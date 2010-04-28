#!/usr/bin/env python
'''
This script is used to create a backbone project
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of project.
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
