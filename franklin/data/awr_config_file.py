## assemble to reference configuration file
# This is a python dictionary that contains three keys:
# .- working directory that is mandatory
# .- reference that is a python dictionary. Mandatory!!
# .- reads. A list of pythond dictionary

{'working_dir' : 'working_dir',
 'reference'   : {'name'      : 'reference',
                   'format'    : 'fasta',
                   'seq_fname' : 'reference.fasta',
                   'qual_fname': 'refernce.qual.fasta'},

 'reads'       : [{'name'       :'solexa1',
                   'seq_fname'  : 'name1',
                   'qual_fname' : 'name1',
                   'format'     : 'format1',
                   'seq_type'   : 'solexa',
                   'aligner'    : 'nucmer',
                   'aligner_parameters' : '-l 20 -c 20'},

                  {'name'       :'sanger1',
                   'seq_fname'  : 'name2',
                   'qual_fname' : 'name1',
                   'format'     : 'fasta',
                   'seq_type'   : 'sanger',
                   'aligner'    : 'nucmer',
                   'aligner_parameters' : ''}
                 ]
}