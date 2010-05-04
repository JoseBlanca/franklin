'''
Created on 12/03/2010

@author: peio
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# franklin is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.


from franklin.seq.seq_annotation import (create_cdna_intron_annotator,
                                         create_ortholog_annotator,
                                         create_description_annotator,
                                         create_microsatellite_annotator,
                                         create_orf_annotator,
                                         create_go_annotator)

annotate_cdna_introns = {'function': create_cdna_intron_annotator,
                         'arguments':{'genomic_db':None,
                                      'genomic_seqs_fhand':None},
                         'type':'mapper',
                         'name':'annotate_cdna_introns',
            'comment': 'It annotates introns comparing with a reference genome'}

annotate_orthologs = {'function': create_ortholog_annotator,
                      'arguments':{'blast':None, 'reverse_blast':None,
                                   'species': None},
                      'type':'mapper' ,
                      'name':'annotate_orthologs',
                      'comment': 'It annotates orthologs using reverse blasts'}

annotate_with_descriptions = {'function': create_description_annotator,
                              'arguments':{'blasts':None},
                              'type':'mapper' ,
                              'name':'annotate_descriptions',
           'comment': 'It annotates using the first hit of the given databases'}

annotate_microsatellites = {'function': create_microsatellite_annotator,
                        'arguments':{},
                        'type':'mapper' ,
                        'name':'annotate_microsatellites',
                        'comment': 'It annotates The microsatellites'}

annotate_orfs = {'function': create_orf_annotator,
                 'arguments':{'parameters':None},
                 'type':'mapper' ,
                 'name':'annotate_orfs',
                 'comment': 'It annotates The orf'}

annotate_gos = {'function': create_go_annotator,
                 'arguments':{'blast':None},
                 'type':'mapper' ,
                 'name':'annotate_gos',
                 'comment': 'It annotates the gos'}


ANNOT_STEPS = [annotate_cdna_introns, annotate_orthologs,
               annotate_with_descriptions, annotate_microsatellites,
               annotate_orfs, annotate_gos]
