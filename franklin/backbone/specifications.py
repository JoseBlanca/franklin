'''
Created on 12/03/2010

@author: peio
'''
BACKBONE_DIRECTORIES = {
    'config_file': 'franklin.conf',
    'external_software_config': 'config_data',
    'original_reads': 'reads/original',
    'cleaned_reads': 'reads/cleaned',
    'assembly_input': 'assembly/input',
    'assemblies': ('assembly', ''),
    'assembly_result': ('assembly', 'result'),
    'mappings': ('mapping', ''),
    'mapping_result': ('mapping', 'result'),
    'mapping_reference': 'mapping/reference',
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
    'error_dir': 'franklin_errors',
    'go_files': 'annotations/go_files',
                       }
BACKBONE_BASENAMES = {
    'contigs':'contigs',
    'mapping_reference':'reference',
    'merged_bam':'merged.bam',
    'snv_result':'all.snvs',
    'merged_frg':'all_seq.frg',
    'blast_basename':'blast'
}
