'''
Created on 2009 api 7

@author: peio
'''
from biolib.SeqVariation import seqvariations_in_alignment
from biolib.cafparser import CafParser
#from sys import argv 

def main():
    ''' Main function where we find snps'''
    #caf_file = argv[1]
    #caf_file = "/home/peio/work_in/assemblers/eucalyptus_c5476.caf"
    caf_file = "/home/peio/assemblers/mira_exec/eucalyptus_d_results/eucalyptus_out.caf"
    #caf_file = '/home/jose/tmp/example.caf'
    caf_parser = CafParser(caf_file)
    var_count = {}
    contig_count = 0
    for contig in caf_parser.contigs():
        contig_count += 1
        for seqvar in seqvariations_in_alignment(contig):
            kind = seqvar.kind()
            if not kind in var_count:
                var_count[kind] = 0
            var_count[kind] += 1
    print 'contigs: ' + str(contig_count)
    print 'variations:' + str(var_count)


if __name__ == '__main__':
    main()
