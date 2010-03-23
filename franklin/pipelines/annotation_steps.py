'''
Created on 12/03/2010

@author: peio
'''
from franklin.seq.seq_annotation import (create_cdna_intron_annotator,
                                         create_ortholog_annotator,
                                         create_description_annotator)

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

annotate_with_descriptions = {'function': create_description_annotator,
                        'arguments':{'blasts':None},
                        'type':'mapper' ,
                        'name':'annnoatate_descriptions',
           'comment': 'It annotates using the first hit of the given databases'}
