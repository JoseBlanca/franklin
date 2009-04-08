'''
Created on 2009 api 7

@author: peio
'''
from biolib.SeqVariation import seqvariations_in_alignment
from biolib.cafparser import CafParser
from sys import argv 

def main():
    ''' Main function where we find snps'''
    #caf_file = argv[1]
    caf_file = "/home/peio/work_in/assemblers/eucalyptus_c5476.caf"
    print "Creating the index"
    caf_parser = CafParser(caf_file)
    print "Index Finished"
    for contig in caf_parser.contigs():
#        print contig
        for seqvar in seqvariations_in_alignment(contig):
            print "Seqvar is", seqvar.kind() 
            print seqvar
            


if __name__ == '__main__':
    main()