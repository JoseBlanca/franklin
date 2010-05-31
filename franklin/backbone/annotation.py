'''
This module is part the ngs_backbone. It performas analyses related to sequence
annotation.

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

import shutil, os, copy
from os.path import join, exists
from franklin.backbone.analysis import Analyzer
from tempfile import NamedTemporaryFile
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.backbone.blast_runner import (backbone_blast_runner,
                                            make_backbone_blast_db,
                                            guess_blastdb_kind)
from franklin.pipelines.annotation_steps import (annotate_cdna_introns,
                                                 annotate_orthologs,
                                                 annotate_with_descriptions,
                                                 annotate_microsatellites,
                                                 annotate_orfs,
                                                 annotate_gos)
from franklin.sam import create_bam_index
from franklin.pipelines.snv_pipeline_steps import (
                                            unique_contiguous_region_filter,
                                            close_to_intron_filter,
                                            high_variable_region_filter,
                                            close_to_snv_filter,
                                            close_to_limit_filter,
                                            major_allele_freq_filter,
                                            kind_filter,
                                            cap_enzyme_filter,
                                            is_variable_filter,
                                            ref_not_in_list)
from franklin.seq.writers import (SequenceWriter, GffWriter, SsrWriter,
                                  OrfWriter, OrthologWriter)
from franklin.snv.writers import VariantCallFormatWriter, SnvIlluminaWriter

from franklin.utils.misc_utils import VersionedPath, xml_itemize
from franklin.utils.cmd_utils import b2gpipe_runner

class AnnotationAnalyzer(Analyzer):
    'It annotates the introns in cdna sequences'

    def _get_inputs_and_prepare_outputs(self):
        'It creates the output dir and it returns the inputs'
        output_dirs = self._create_output_dirs()
        inputs = self._get_input_fpaths()
        return inputs, output_dirs

    def _run_annotation(self, pipeline, configuration, inputs, output_dir):
        'It runs the analysis.'


        self._log({'analysis_started':True})
        repr_fpaths = inputs['repr']
        try:
            seqs_fpaths = inputs['input']
        except KeyError:
            seqs_fpaths = []
        seqs_paths = self._get_seq_or_repr_path(seqs_fpaths, repr_fpaths)
        for seq_path in seqs_paths:
            seq_fpath = seq_path.last_version
            temp_repr = NamedTemporaryFile(suffix='.repr', mode='a',
                                           delete=False)
            in_fhands = {'in_seq': open(seq_fpath)}

            writer = SequenceWriter(fhand=temp_repr,
                                    file_format='repr')

            if seq_path.basename in configuration:
                #there is a different configuration for every file to annotate
                config = configuration[seq_path.basename]
            else:
                config = configuration

            seq_pipeline_runner(pipeline, configuration=config,
                                in_fhands=in_fhands,
                                processes=self.threads,
                                writers={'repr': writer})
            temp_repr.close()
            repr_path = VersionedPath(os.path.join(output_dir,
                                                   seq_path.basename + '.repr'))
            repr_fpath = repr_path.next_version
            shutil.move(temp_repr.name, repr_fpath)
        self._log({'analysis_finished':True})

    def _get_seq_or_repr_path(self, seqs_paths, repr_paths):
        'It returns for every file the repr or the seq file'
        repr_paths_ = dict([(path.basename, path) for path in repr_paths])
        if not seqs_paths:
            return repr_paths
        new_seq_paths = []
        for path in seqs_paths:
            basename = path.basename
            if basename in repr_paths_:
                new_seq_paths.append(repr_paths_[basename])
            else:
                new_seq_paths.append(path)
        return new_seq_paths

class AnnotateOrthologsAnalyzer(AnnotationAnalyzer):
    '''It annotates the orthologs for the working species against the species
    given in the configuration file'''

    def run(self):
        'It runs the analysis'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        output_dir = output_dirs['result']
        blast_settings = self._project_settings['blast']
        settings = self._project_settings['Annotation']['ortholog_annotation']
        ortholog_databases = settings['ortholog_databases']


        general_settings = self._project_settings['General_settings']
        project_dir = general_settings['project_path']

        #first we need some blasts
        blasts = {}
        for input_ in inputs['input']:
            for database in ortholog_databases:
                if 'kind' in blast_settings[database]:
                    db_kind = blast_settings[database]['kind']
                else:
                    db_kind = guess_blastdb_kind(blast_settings[database]['path'])

                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'blastx'
                blastdb = blast_settings[database]['path']
                #this could be different adding something to the settings
                blastdb_seq_fpath = blastdb
                blast = backbone_blast_runner(query_fpath=input_.last_version,
                                              project_dir=project_dir,
                                              blast_program=blast_program,
                                              blast_db=blastdb,
                                              dbtype=db_kind,
                                              threads=self.threads)
                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'tblastn'
                reverse_blast = backbone_blast_runner(
                                              query_fpath=blastdb_seq_fpath,
                                              project_dir=project_dir,
                                              blast_program=blast_program,
                                              blast_db_seq=input_.last_version,
                                              dbtype='nucl',
                                              threads=self.threads)
                if input_ not in blasts:
                    blasts[input_] = {}
                blasts[input_][database] = {'blast':blast,
                                            'reverse_blast':reverse_blast}

        pipeline = []
        configuration = {}
        for database in ortholog_databases:
            step = copy.deepcopy(annotate_orthologs)
            step['name_in_config'] = database
            #an annotation step for every ortholog database
            pipeline.append(step)
            for input_ in inputs['input']:
                reverse_blast = ''
                step_config = {
                    'blast':{'blast': blasts[input_][database]['blast']},
                    'reverse_blast':{'blast':
                                     blasts[input_][database]['reverse_blast']},
                    'species': database}
                if input_.basename not in configuration:
                    configuration[input_.basename] = {}
                configuration[input_.basename][database] = step_config

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

class AnnotateIntronsAnalyzer(AnnotationAnalyzer):
    'It annotates the introns in cdna sequences'

    def run(self):
        'It runs the analysis.'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        repr_dir = output_dirs['repr_dir']
        pipeline = [annotate_cdna_introns]

        general_settings = self._project_settings['General_settings']

        settings = self._project_settings['Annotation']
        genomic_seqs = settings['Cdna_intron_annotation']['genomic_seq_file']
        if settings['Cdna_intron_annotation']['genomic_db'] is None:
            project_path = general_settings['project_path']
            genomic_db = make_backbone_blast_db(project_path, genomic_seqs,
                                                dbtype='nucl')
        else:
            genomic_db = settings['Cdna_intron_annotation']['genomic_db']
        configuration = {'annotate_cdna_introns': {'genomic_db':genomic_db,
                                             'genomic_seqs_fhand':genomic_seqs}}

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=repr_dir)

class SnvCallerAnalyzer(AnnotationAnalyzer):
    'It performs the calling of the snvs in a bam file'

    def run(self):
        'It runs the analysis.'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        output_dir = output_dirs['result']
        merged_bam = inputs['merged_bam']
        create_bam_index(merged_bam.last_version)

        pipeline = 'snv_bam_annotator'
        bam_fpath = merged_bam.last_version
        configuration = {'snv_bam_annotator': {'bam_fhand':bam_fpath}}
        settings = self._project_settings
        if 'Snvs' in settings:
            snv_settings = settings['Snvs']
            # read egde conf
            read_edge_conf = self._configure_read_edge_conf(snv_settings)
            configuration['snv_bam_annotator']['read_edge_conf'] = read_edge_conf
            for config_param in ('min_quality', 'min_mapq', 'min_num_alleles'):
                param_value = int(snv_settings[config_param])
                configuration['snv_bam_annotator'][config_param] = param_value
        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=output_dir)

    @staticmethod
    def _configure_read_edge_conf(snv_settings):
        '''It takes from the settings the edges to take into account in snv
        caller'''
        read_edge_conf = {}
        if 'edge_removal' in snv_settings:
            for pl_side, edge in snv_settings['edge_removal'].items():
                platform, side = pl_side.split('_')
                if platform not in read_edge_conf:
                    read_edge_conf[platform] = [None, None]
                side = 0 if side == 'left' else 1
                read_edge_conf[platform][side] = edge
        return read_edge_conf

class WriteAnnotationAnalyzer(Analyzer):
    'It writes all the output annotation files'
    def run(self):
        'It runs the analysis.'
        output_dir = self._create_output_dirs()['result']
        inputs = self._get_input_fpaths()
        repr_paths = inputs['repr']


        output_files = {'vcf': ('vcf',),
                        'orf':('orf_seq.fasta', 'orf_pep.fasta'),
                        'ssr':('ssr',),
                        'gff':('gff3',),
                        'orthologs':('orthologs',),}

        for seq_path in repr_paths:
            outputs = {}
            for kind, extensions in output_files.items():
                outputs[kind] = []
                for extension in extensions:
                    output_fpath = os.path.join(output_dir,
                                            seq_path.basename + '.' + extension)
                    if os.path.exists(output_fpath):
                        os.remove(output_fpath)
                    output_fhand = open(output_fpath, 'a')
                    outputs[kind].append(output_fhand)

            for kind, output in outputs.items():
                if len(output) == 1:
                    outputs[kind] = output[0]

            in_fhands = {'in_seq': open(seq_path.last_version)}

            writers = {}
            if 'repr' in outputs:
                writers['repr'] = SequenceWriter(fhand=outputs['repr'],
                                                 file_format='repr')

            if 'vcf' in outputs:
                ref_name = os.path.basename(in_fhands['in_seq'].name)
                fhand = outputs['vcf']
                grouping = self._project_settings['Snvs']['vcf_grouping']
                writers['vcf'] = VariantCallFormatWriter(fhand=fhand,
                                                        reference_name=ref_name,
                                                        grouping=grouping)
            if 'gff' in outputs:
                default_type = None
                writers['gff'] = GffWriter(fhand=outputs['gff'],
                                           default_type=default_type)
            if 'orf' in outputs:
                fhand, pep_fhand = outputs['orf']
                writers['orf'] = OrfWriter(fhand=fhand, pep_fhand=pep_fhand)

            if 'ssr' in outputs:
                writers['ssr'] = SsrWriter(fhand=outputs['ssr'])

            if 'orthologs' in outputs:
                writers['orthologs'] = \
                                      OrthologWriter(fhand=outputs['orthologs'])

            if 'snv_illumina' in outputs:
                writers['snv_illumina'] = \
                                SnvIlluminaWriter(fhand=outputs['snv_illumina'])

            feature_counter = seq_pipeline_runner(pipeline=None,
                                                  configuration=None,
                                                  in_fhands=in_fhands,
                                                  writers = writers)

            # We need to close fhands and remove void files.
            # sequence writer could have a qual fhand
            # orf writer has a pep_fhand
            for kind, fhands in outputs.items():
                kind_key = 'sequence' if kind == 'quality' else kind
                if kind in feature_counter:
                    self._close_and_remove_files(fhands,
                                                 feature_counter[kind_key])

    @staticmethod
    def _close_and_remove_files(fhands, num_features):
        'it closes and removes files that are empty'
        if not isinstance(fhands, tuple) and not isinstance(fhands, list):
            fhands = (fhands,)
        for fhand in fhands:
            fpath = fhand.name
            fhand.close()
            if not num_features and os.path.exists(fpath):
                os.remove(fpath)

class SnvFilterAnalyzer(AnnotationAnalyzer):
    'It performs the filtering analysis of the snvs'

    def run(self):
        'It runs the analysis.'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        repr_dir = output_dirs['repr_dir']
        pipeline, configuration = self._get_snv_filter_specs()

        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=repr_dir)
    def _get_snv_filter_specs(self):
        'It gets the pipeline and configuration from settings'
        configuration = {}
        settings = self._project_settings['snv_filters']
        for filter_data in settings.values():
            if not filter_data['use']:
                continue
            name = filter_data['name']
            unique_name = filter_data['unique_name'] if 'unique_name' in filter_data else name
            filter_config = {}
            for argument, value in filter_data.items():
                if argument == 'use' or argument == 'name':
                    continue
                if name == 'ref_not_in_list':
                    filter_config['seq_list'] = []
                    for item  in open(value):
                        filter_config['seq_list'].append(item.strip())
                else:
                    filter_config[argument] = value
                if name == 'by_kind' and argument == 'kind':
                    filter_config[argument] =  self._snv_kind_to_franklin(value)

            filter_config['name'] = name
            configuration[unique_name] = filter_config
        pipeline = self._get_pipeline_from_step_name(configuration)

        return pipeline, configuration

    @staticmethod
    def _snv_kind_to_franklin(value):
        '''it translates the snv kind value to the number code understanded by
        franklin'''
        if value in (0, '0', 'SNP', 'snp'):
            return 0
        elif value in (3, '3', 'INVARIANT', 'invariant'):
            return 3
        elif value in (4, '4', 'INDEL', 'indel'):
            return 4

    def _get_pipeline_from_step_name(self, configuration):
        'It gets the pipeline steps from filter_names'
        translator = {'uniq_contiguous' : unique_contiguous_region_filter,
                      'close_to_intron':close_to_intron_filter,
                      'high_variable_region': high_variable_region_filter,
                      'close_to_snv':close_to_snv_filter,
                      'close_to_limit': close_to_limit_filter,
                      'maf': major_allele_freq_filter,
                      'by_kind' : kind_filter,
                      'cap_enzyme': cap_enzyme_filter,
                      'is_variable': is_variable_filter,
                      'ref_not_in_list':ref_not_in_list}
        pipeline = []
        for unique_name in configuration:
            step_name = configuration[unique_name]['name']
            #we have to remove the step name from the configuration because this
            #dict will be used as **kargs
            if 'unique_name' in configuration[unique_name]:
                del configuration[unique_name]['unique_name']
            del configuration[unique_name]['name']
            step = copy.deepcopy(translator[step_name])
            step['name_in_config'] = unique_name
            pipeline.append(step)
        return pipeline

class AnnotateDescriptionAnalyzer(AnnotationAnalyzer):
    '''This class is used to annotate sequences using the first hit on the given
    blast databases'''
    def run(self):
        'It runs the analysis'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        repr_dir = output_dirs['repr_dir']
        blast_settings = self._project_settings['blast']

        settings = self._project_settings['Annotation']
        annot_settings = settings['description_annotation']
        description_databases = annot_settings['description_databases']

        general_settings = self._project_settings['General_settings']

        #first we need some blasts
        project_dir = general_settings['project_path']
        blasts = {}
        for input_ in inputs['input']:
            input_fpath = input_.last_version
            for database in description_databases:
                if 'kind' in blast_settings[database]:
                    db_kind = blast_settings[database]['kind']
                else:
                    db_kind = guess_blastdb_kind(blast_settings[database]['path'])

                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'blastx'

                blastdb = blast_settings[database]['path']
                blast = backbone_blast_runner(query_fpath=input_fpath,
                                                project_dir=project_dir,
                                                blast_program=blast_program,
                                                blast_db=blastdb,
                                                dbtype=db_kind,
                                                threads=self.threads)
                if input_ not in blasts:
                    blasts[input_fpath] = []
                blasts[input_fpath].append({'blast':blast, 'modifier':None})
        #print blasts
        pipeline = []
        configuration = {}
        for database in description_databases:
            step = annotate_with_descriptions
            step['name_in_config'] = database
            pipeline.append(step)
            for input_ in inputs['input']:
                step_config = {'blasts': blasts[input_.last_version]}
                configuration[input_.basename] = {}
                configuration[input_.basename][database] = step_config
        #print configuration
        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=repr_dir)

class AnnotateMicrosatelliteAnalyzer(AnnotationAnalyzer):
    '''This class is used to annotate the microsatellite of a sequence as a
    feature'''

    def run(self):
        'It runs the analysis.'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        repr_dir = output_dirs['repr_dir']

        pipeline = [annotate_microsatellites]
        configuration = {}
        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=repr_dir)

class AnnotateOrfAnalyzer(AnnotationAnalyzer):
    '''This class is used to annotate the orf of a sequence as a feature'''

    def run(self):
        'It runs the analysis.'
        matrix = self._project_settings['Annotation']['orf_annotation']['estscan_matrix']
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        repr_dir = output_dirs['repr_dir']
        pipeline = [annotate_orfs]

        configuration = {'annotate_orfs': {'parameters':{'matrix':matrix}}}
        return self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=repr_dir)

class AnnotateGoAnalyzer(AnnotationAnalyzer):
    'This class is used to annotate gos using blast2go'

    def run(self):
        'It runs the analysis.'
        inputs, output_dirs = self._get_inputs_and_prepare_outputs()
        repr_dir            = output_dirs['repr_dir']
        result_dir          = output_dirs['result']
        blast_settings      = self._project_settings['blast']

        go_settings    = self._project_settings['Annotation']['go_annotation']

        go_database = go_settings['blast_database']
        create_dat  = go_settings['create_dat_file']
        java_memory = go_settings['java_memory']
        prop_fpath  = go_settings['b2g_properties_file']

        blast2go = {}
        for input_ in inputs['input']:
            input_fpath = input_.last_version
            annot_fpath = join(result_dir, input_.basename + '.b2g.annot')
            if not exists(annot_fpath):
                go_blast_settings = blast_settings[go_database]
                blast = self._get_b2g_blast(input_fpath, go_blast_settings)
                if create_dat:
                    dat_fpath = join(result_dir, input_.basename + '.b2g.dat')
                else:
                    dat_fpath = None

                b2gpipe_runner(blast, annot_fpath=annot_fpath,
                               dat_fpath=dat_fpath, java_memory=java_memory,
                               prop_fpath=prop_fpath)
                blast.close()

            blast2go[input_fpath] = annot_fpath

        # prepare pipeline
        pipeline = [annotate_gos]
        configuration = {}

        for input_ in  inputs['input']:
            input_fpath = input_.last_version

            step_config = {'annot_fpath': blast2go[input_fpath]}

            configuration[input_.basename] = {'annotate_gos': step_config}
        result = self._run_annotation(pipeline=pipeline,
                                      configuration=configuration,
                                      inputs=inputs, output_dir=repr_dir)

        return result

    def _get_b2g_blast(self, input_fpath, goblast_settings):
        'It gets a chopped blast ready for use with blast2go'
        if 'kind' in goblast_settings:
            db_kind = goblast_settings['kind']
        else:
            db_kind = guess_blastdb_kind(goblast_settings['path'])

        if db_kind == 'nucl':
            blast_program = 'tblastx'
        else:
            blast_program = 'blastx'

        blastdb     = goblast_settings['path']

        project_dir = self._project_settings['General_settings']['project_path']
        blast       = backbone_blast_runner(query_fpath=input_fpath,
                                            project_dir=project_dir,
                                            blast_program=blast_program,
                                            blast_db=blastdb,
                                            dbtype=db_kind,
                                            threads=self.threads)

        chop_big_xml, num_items = True, 2
        if chop_big_xml:
            #chopped_blast = open('/tmp/blast_itemized.xml', 'w')
            chopped_blast = NamedTemporaryFile(suffix='.xml')
            for blast_parts in xml_itemize(blast, 'Iteration', num_items):
                chopped_blast.write(blast_parts)
            chopped_blast.flush()
            return chopped_blast
        else:
            return open(blast)

DEFINITIONS = {
    'filter_snvs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'}},
         'analyzer': SnvFilterAnalyzer},
    'annotate_snvs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'merged_bam':
                {'directory':'mapping_result',
                 'file': 'merged_bam'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': SnvCallerAnalyzer},
    'annotate_introns':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateIntronsAnalyzer},
    'write_annotations':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_result'}},
         'analyzer': WriteAnnotationAnalyzer},
    'annotate_orthologs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'result':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateOrthologsAnalyzer},
    'annotate_descriptions':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateDescriptionAnalyzer},

    'annotate_microsatellites':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateMicrosatelliteAnalyzer},
    'annotate_orfs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'}},
         'analyzer': AnnotateOrfAnalyzer},
    'annotate_gos':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            'input':
                {'directory': 'annotation_input',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'},
                    'result':{'directory':'annotation_result'}},
         'analyzer': AnnotateGoAnalyzer},
    }
