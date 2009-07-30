#!/usr/bin/env python
'''
This script search snps in an assembly. The assembly could be in ace or caf
format

Created on 2009 mai 4
@author: peio
'''
from biolib.seqvariation import seqvariations_in_alignment, seqvar_summary
from biolib.contig_io import get_parser_by_name
from biolib.contig_cleaner import contig_strip, water_alignment_strip
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

    if options.soutfile is None:
        soutfile = "".join(options.infile.split('.')[:-1]) + '.snp_info'
        io_fhands['snp_outfile'] = open(basename(soutfile), 'w')
    else:
        io_fhands['snp_outfile'] = open(options.soutfile, 'w')

    if options.coutfile is None:
        coutfile = "".join(options.infile.split('.')[:-1]) + '.contig_info'
        io_fhands['contig_outfile'] = open(basename(coutfile), 'w')
    else:
        io_fhands['contig_outfile'] = open(options.coutfile, 'w')
    if options.c_strip is None:
        c_strip = False
    else:
        c_strip = True

    return io_fhands, c_strip

def main():
    ''' Main function where we find snps'''
    io_fhands, c_strip = set_parameters()
    infhand            = io_fhands['infile']
    snp_fhand          = io_fhands['snp_outfile']
    contig_fhand       = io_fhands['contig_outfile']

    # Get the porper parser for input file format
    parser = get_parser_by_name(infhand.name)

    # write ouput file headers
    snp_fhand.write('#Id\tLoc_start\tLoc_end\tType\tpic\tcap\n')
    contig_fhand.write('format-version:1\n')

    for contig in parser.contigs():
        contig_name = contig.consensus.sequence.name
        if c_strip:
            contig = contig_strip(contig, 13) # este numero es aproximado
        contig = water_alignment_strip(contig)
        print ".-Searching in contig: %s" % contig_name
        var_count = 0

        for seqvar in seqvariations_in_alignment(contig):
            snp_print = seqvar_summary(seqvar, contig_name)
            snp_fhand.write(snp_print)

            var_count += 1
        if var_count == 0:
            percentaje_snp = 0
        else:
            percentaje_snp = (var_count/float(len(contig.consensus))) * 100

        toprint = '%s\t%d\t%f\n' % (contig_name, var_count, percentaje_snp)
        contig_fhand.write(toprint)

    print "Finished"


if __name__ == '__main__':
    main()
