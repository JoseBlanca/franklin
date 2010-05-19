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

BACKBONE_DIRECTORIES = {
    'config_file': 'backbone.conf',
    'log_file': 'backbone.log',
    'external_software_config': 'config_data',
    'original_reads': 'reads/original',
    'cleaned_reads': 'reads/cleaned',
    'assembly_input': 'assembly/input',
    'assemblies': ('assembly', ''),
    'assembly_result': ('assembly', 'result'),
    'mappings': ('mapping', ''),
    'mapping_result': ('mapping', 'result'),
    'mapping_reference': 'mapping/reference',
    'mapping_stats':('mapping', 'result/stats'),
    'mappings_by_readgroup': ('mapping', 'result/by_readgroup'),
    'pileups':'mapping/result/pileups',
    'snvs':'annotations/snvs',
    'info':'info',
    'original_reads_stats': 'reads/original/stats',
    'cleaned_reads_stats': 'reads/cleaned/stats',
    'annotation_repr':'annotations/repr',
    'annotation_input':'annotations/input',
    'annotation_result':'annotations/result',
    'blast_dir':'annotations/blast',
    'blast_databases':'annotations/blast/databases',
    'error_dir': 'backbone_errors',
                       }
BACKBONE_BASENAMES = {
    'contigs':'contigs',
    'mapping_reference':'reference',
    'merged_bam':'merged.bam',
    'snv_result':'all.snvs',
    'merged_frg':'all_seq.frg',
    'blast_basename':'blast'
}
