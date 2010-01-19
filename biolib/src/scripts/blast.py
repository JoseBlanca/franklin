'''
Created on 19/01/2010


@author: peio
'''
from optparse import OptionParser

''def parse_options():
    'It parses the command line arguments'
    parser = OptionParser()
    parser.add_option('-d', dest='database', help='database')
    parser.add_option('-p', dest='program', help='blast program')
    parser.add_option('-e', dest="expect", help='evalue')
    parser.add_option('-v', dest='nhitsv', help='nhitsv')
    parser.add_option('-b', dest='nhitsb', help='nhitsb')
    parser.add_option('-m', dest="outformat", help='output file format')
    parser.add_option('-i', dest='query', help='Input file')

    return parser

def set_parameters():
    'Set parameters'
    # Set parameters
    parser  = parse_options()
    options = parser.parse_args()[0]

    database = options.database
    program  = options.program
    expect   = options.expect
    nhitsv   = options.nhitsv
    nhitsb   = options.nhitsb
    outformat = options.outformat
    query    = options.query
    return database, program, expect, nhitsv, nhitsb, outformat, query

def blast():
    "lucy runner"
    database, program, expect, nhitsv, nhitsb, outformat, query = \
                                                                set_parameters()
    cmd = [program, '-db', database, '-query', query, '-outfmt', outformat]
    if expect is not None:
        cmd.extend(['-evalue', expect])
    if nhitsv is not None:
        cmd.extend(['-num_descriptions', nhitsv])
    if nhitsb is not None:
        cmd.extend(['-num_alignments', nhitsb])

    stdout = call(cmd, raise_on_error=True)[0]
    print stdout







if __name__ == '__main__':
   blast()
