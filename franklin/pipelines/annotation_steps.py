'''
Created on 12/03/2010

@author: peio
'''
from franklin.seq.seq_annotation import (create_cdna_intron_annotator,
                                         create_ortholog_annotator,)

annotate_cdna_introns = {'function': create_cdna_intron_annotator,
                        'arguments':{'genomic_db':None,
                                     'genomic_seqs_fhand':None},
                        'type':'mapper' ,
                        'name':'annnoatate_cdna_introns',
            'comment': 'It annotates introns comparing with a reference genome'}

annotate_orthologs = {'function': create_ortholog_annotator,
                       'arguments':{'blast':None, 'reverse_blast':None,
                                    'species': None},
                       'type':'mapper' ,
                       'name':'annnoatate_orthologs',
                       'comment': 'It annotates orthologs using reverse blasts'}
