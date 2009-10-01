#!/usr/bin/env python
'''
This script search snps in an assembly. The assembly could be in ace or caf
format

Created on 2009 mai 4
@author: peio
'''
from biolib.seqvariation import (seqvars_in_contigs, seqvar_summary)
from biolib.contig_io import get_parser
from biolib.pipelines import pipeline_runner

from optparse import OptionParser
from os.path import basename

def parse_options():
    'It parses the command line arguments'
    parser = OptionParser('usage: %prog -i infile [-t]...', version='%prog 1.0')
    parser.add_option('-s', '--soutfile',  dest='soutfile',
                      help='snp information output file')
    parser.add_option('-c', '--coutfile',  dest='coutfile',
                      help='contig information output file')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Input file')
    parser.add_option('-t', '--contig_strip', dest = 'c_strip',
                     action='store_true', help = 'Strip contig?'  )
    return  parser.parse_args()

def set_parameters():
    '''It sets the parameters for the script.

    For some options there are some default values that will be applied here.'''
    options = parse_options()[0]
    io_fhands = {}

    if options.infile is None:
        raise RuntimeError('Script at least needs an input file (caf|ace)')
    else:
        io_fhands['infile'] = open(options.infile, 'r')

#    if options.soutfile is None:
#        soutfile = "".join(options.infile.split('.')[:-1]) + '.snp_info'
#        io_fhands['snp_outfile'] = open(basename(soutfile), 'w')
#    else:
#        io_fhands['snp_outfile'] = open(options.soutfile, 'w')
#
#    if options.coutfile is None:
#        coutfile = "".join(options.infile.split('.')[:-1]) + '.contig_info'
#        io_fhands['contig_outfile'] = open(basename(coutfile), 'w')
#    else:
#        io_fhands['contig_outfile'] = open(options.coutfile, 'w')

    return io_fhands

def main():
    ''' Main function where we find snps'''
    io_fhands    = set_parameters()
    infhand      = io_fhands['infile']
#    snp_fhand    = io_fhands['snp_outfile']
#    contig_fhand = io_fhands['contig_outfile']

    # Get the porper parser for input file format
    parser = get_parser(infhand)

    #clean and filter non valid contigs/read regions
    contigs = parser.contigs()
    contigs = pipeline_runner('contig_clean', contigs)

    # perform the search and snp filtering
    seq_var_iter = seqvars_in_contigs(contigs)
    seq_var_iter = pipeline_runner('snp_clean', seq_var_iter)

    for seq_var in seq_var_iter:
        print seqvar_summary(seq_var)



if __name__ == '__main__':
    main()
