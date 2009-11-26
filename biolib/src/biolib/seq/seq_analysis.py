'''
Created on 26/11/2009

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of biolib.
# biolib is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# biolib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with biolib. If not, see <http://www.gnu.org/licenses/>.

from biolib.utils.cmd_utils import create_runner
from biolib.alignment_search_result import (FilteredAlignmentResults,
                                            get_alignment_parser)

def _infer_introns_from_match_parts(alignments):
    'Given a match with several match parts it returns the introns'
    alignment = alignments.next()
    match = alignment['matches'][0]
    from pprint import pprint
    pprint(match)


def infer_introns_for_cdna(sequence, genomic_db, blast_db_path):
    '''Doing a blast with the sequences against the genomic db it infers the
    positions of introns'''

    #first we run the blast
    parameters = {'database': genomic_db, 'program':'tblastx'}

    filters = [{'kind'          : 'min_length',
                'min_length_bp' : 15}]
    blast_runner = create_runner(kind='blast', parameters=parameters,
                                 environment={'BLASTDB':blast_db_path})
    blast_fhand = blast_runner(sequence)[0]

    #now we parse the blast
    blast_parser = get_alignment_parser('blast')
    blast_result = blast_parser(blast_fhand)

    # We filter the results with appropiate  filters

    alignments = FilteredAlignmentResults(match_filters=filters,
                                          results=blast_result)

    #now we have to guess the introns
    introns = _infer_introns_from_match_parts(alignments)
    return introns
