#!/usr/bin/env python
'''
It creates a marker centralized file to be used for the rest of the scripts to
normalize all marker name and type fields

@author: peio
'''
from optparse import OptionParser
import os
from franklin.gmod.markers import (write_markers, parse_markersfile,
                                   add_markers_gff3)


def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-i', '--infiles', dest='infiles', help='infiles')
    parser.add_option('-m', '--markersfile', dest='markersfile',
                      help='markers file')
    parser.add_option('-n', '--name_correlations', dest='name_cor',
                      help='Evaluable file with name correlations')
    parser.add_option('-t', '--type_correlations', dest='type_cor',
                      help='Evaluable file with type correlations')
    parser.add_option('-s', '--sequence_correlations', dest='sequence_cor',
                      help='Evaluable file with sequence correlations')
    return parser

def set_parameters():
    'Set parameters'
    parser  = parse_options()
    options = parser.parse_args()[0]

    if options.infiles is None:
        parser.error('In file required')
    else:
        infhands = [open(infile)for infile in options.infiles.split(',')]

    if options.markersfile is None:
        parser.error('We need a markers file to read or write markers to')
    else:
        markersfname =  options.markersfile

    if options.name_cor is None:
        name_cor = None
    else:
        name_cor = correlations_in_file(open(options.name_cor))

    if options.type_cor is None:
        type_cor = None
    else:
        type_cor = correlations_in_file(open(options.type_cor))
    if options.sequence_cor is None:
        sequence_cor = None
    else:
        sequence_cor = correlations_in_file(open(options.sequence_cor))

    return infhands, markersfname, name_cor, type_cor, sequence_cor

def main():
    'The main part'
    # set parameters
    infhands, markersfname, name_cor, type_cor, sequence_cor = set_parameters()

    #read markes file if exists
    if os.path.exists(markersfname):
        markersfhand = open(markersfname)
        markers = parse_markersfile(markersfhand)
        markersfhand.close()
    else:
        markers = {}

    # add files to marker structure
    for infhand in infhands:
        add_markers(infhand, markers, name_cor, type_cor, sequence_cor)

    # write markers to file:
    markersfhand = open(markersfname, 'w')
    write_markers(markersfhand, markers)
    markersfhand.close()

def correlations_in_file(fhand):
    'It indexes a correlation file'
    correlations = {}
    for line in fhand:
        line = line.strip()
        if not line:
            continue
        name, correlation = line.split()
        correlations[name] = correlation
    return correlations

def add_markers(infhand, markers, name_cor=None, type_cor=None,
                sequence_cor=None):
    'It adds markers from fhand correction the names and adding them as alias'
    # guess file type
    file_type = _guess_file_type(infhand)
    if file_type == 'gff3':
        add_markers_gff3(infhand, markers, name_cor, type_cor, sequence_cor)

def _guess_file_type(infhand):
    'It guesses file type'
    first_line = infhand.readline()
    if first_line.strip() == '##gff-version 3':
        file_type = 'gff3'

    infhand.seek(0)
    return file_type

if __name__ == '__main__':
    main()
