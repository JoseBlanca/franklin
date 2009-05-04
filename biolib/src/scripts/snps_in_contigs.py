'''
Created on 2009 mai 4

@author: peio
'''


from biolib.SeqVariation import seqvariations_in_alignment, calculate_pic, \
                                cap_enzime
from biolib.contig_parser import CafParser, AceParser
from biolib.contig_cleaner import contig_strip, water_alignment_strip
from biolib.biolib_utils import get_start_end
from optparse import OptionParser
from os.path import basename


def main():
    ''' Main function where we find snps'''
    
    parser = OptionParser('usage: %prog [-v] -nNODES ...', version='%prog 1.0')
    parser.add_option('-s', '--soutfile',  dest='soutfile', 
                      help='snp information output file')
    parser.add_option('-c', '--coutfile',  dest='coutfile', 
                      help='contig information output file')
    parser.add_option('-i', '--infile', dest='infile',
                      help='Input file')
    parser.add_option('-t', '--contig_strip', dest = 'c_strip', 
                     action='store_true', help = 'Strip contig?'  )
    (options, args) = parser.parse_args()
  

    if options.infile is None:
        parser.error('Script at least needs an input file (caf|ace)')
    else:
        infile = options.infile
        
    if options.soutfile is None:
        soutfile = "".join(infile.split('.')[:-1]) + '.snp_out'
        soutfile = basename(soutfile)
    else:
        soutfile = options.soutfile
    
    soutfileh = open(soutfile, 'wt')
    soutfileh.write('#Id\tLoc_start\tLoc_end\tType\tpic\tcap\n')
    if options.coutfile is None:
        coutfile = "".join(infile.split('.')[:-1]) + '.contig_out'
        coutfile = basename(coutfile)
    else:
        coutfile = options.coutfile
    coutfileh = open(coutfile, 'wt')
    coutfileh.write('#Id\t%snps in contig\n')
      
    print "Starting file indexing"  
    if infile[-3:].lower() == 'ace':
        print "File type: ace"
        parser = AceParser(infile)
    elif infile[-3:].lower() == 'caf':
        print "File type: caf"
        parser = CafParser(infile)
    print "File indexing finished"
        
    
    for contig in parser.contigs():
        contig_name = contig.consensus.sequence.name
        if options.c_strip is not None:
            contig = contig_strip(contig, 70) # este numero es aproximado
        contig = water_alignment_strip(contig)
        print ".-Searching in contig: %s" % contig_name
        var_count = 0
        for seqvar in seqvariations_in_alignment(contig):
            kind      = seqvar.kind()
            loc_start, loc_end = get_start_end(seqvar.location)
            id_snp    = contig_name + '_' +  str(loc_start)
            pik       = calculate_pic(seqvar)
            enzymes   = cap_enzime(seqvar)
            soutfileh.write("%s\t%d\t%d\t%s\t%f\t%s\n"  % (id_snp, loc_start,
                                                            loc_end, kind, pik, 
                                                            ",".join(enzymes)))
            var_count +=1
        if var_count == 0:
            percentaje_snp = 0
        else:
            percentaje_snp = (var_count/float(len(contig.consensus))) * 100
        coutfileh.write('%s\t%f\n' % (contig_name, percentaje_snp))
        
    print "Finished"    


if __name__ == '__main__':
    main()