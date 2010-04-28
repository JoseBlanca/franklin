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

import shutil, os
from franklin.backbone.analysis import Analyzer
from tempfile import NamedTemporaryFile
from franklin.pipelines.pipelines import seq_pipeline_runner
from franklin.backbone.blast_runner import (backbone_blast_runner,
                                            make_backbone_blast_db)
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
from franklin.utils.misc_utils import VersionedPath, xml_itemize

class AnnotationAnalyzer(Analyzer):
    'It annotates the introns in cdna sequences'

    def _get_inputs_and_prepare_outputs(self):
        'It creates the output dir and it returns the inputs'
        output_dirs   = self._create_output_dirs()
        inputs       = self._get_input_fpaths()
        return inputs, output_dirs

    def _run_annotation(self, pipeline, configuration, inputs, output_dir):
        'It runs the analysis.'
        self._log({'analysis_started':True})
        repr_fpaths  = inputs['repr']
        try:
            seqs_fpaths  = inputs['input']
        except KeyError:
            seqs_fpaths = []
        seqs_paths  = self._get_seq_or_repr_path(seqs_fpaths, repr_fpaths)
        for seq_path in seqs_paths:
            seq_fpath = seq_path.last_version
            temp_repr = NamedTemporaryFile(suffix='.repr', mode='a',
                                           delete=False)
            io_fhands = {'in_seq': open(seq_fpath),
                         'outputs':{'repr':temp_repr}}

            if seq_path.basename in configuration:
                #there is a different configuration for every file to annotate
                config = configuration[seq_path.basename]
            else:
                config = configuration

            seq_pipeline_runner(pipeline, configuration=config,
                                io_fhands=io_fhands)
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

        #first we need some blasts
        project_dir = self._project_settings['General_settings']['project_path']
        blasts = {}
        for input_ in inputs['input']:
            for database in ortholog_databases:
                db_kind = blast_settings[database]['kind']
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
                                              dbtype=db_kind)
                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'tblastn'
                reverse_blast = backbone_blast_runner(
                                              query_fpath=blastdb_seq_fpath,
                                              project_dir=project_dir,
                                              blast_program=blast_program,
                                              blast_db_seq=input_.last_version,
                                              dbtype='nucl')
                if input_ not in blasts:
                    blasts[input_] = {}
                blasts[input_][database] = {'blast':blast,
                                            'reverse_blast':reverse_blast}

        pipeline = []
        configuration = {}
        for database in ortholog_databases:
            step = annotate_orthologs
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
        settings = self._project_settings['Annotation']
        if 'Cdna_intron_annotation' not in settings:
            msg = 'You should set up genomic_db and genomic_seqs in settings'
            raise RuntimeError(msg)

        genomic_seqs = settings['Cdna_intron_annotation']['genomic_seqs']
        if (not 'genomic_db' in settings['Cdna_intron_annotation'] or
            not  settings['Cdna_intron_annotation']['genomic_db']):
            project_path = self._project_settings['General_settings']['project_path']
            genomic_db = make_backbone_blast_db(project_path,
                                                genomic_seqs,
                                                dbtype='nucl')
        else:
            genomic_db = settings['Cdna_intron_annotation']['genomic_db']
        configuration = {'annotate_cdna_introns': {'genomic_db':genomic_db,
                                       'genomic_seqs_fhand':open(genomic_seqs)}}

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
        merged_bam   = inputs['merged_bam']
        create_bam_index(merged_bam.last_version)

        pipeline = 'snv_bam_annotator'
        bam_fhand = open(merged_bam.last_version)
        configuration = {'snv_bam_annotator': {'bam_fhand':bam_fhand}}
        settings = self._project_settings
        if 'Snvs' in settings:
            snv_settings = settings['Snvs']
            # read egde conf
            read_edge_conf = self._configure_read_edge_conf(snv_settings)
            configuration['snv_bam_annotator']['read_edge_conf'] = read_edge_conf
            for confif_param in ('min_quality', 'min_mapq', 'min_num_alleles'):
                if confif_param in snv_settings:
                    param_value = int(snv_settings[confif_param])
                    configuration['snv_bam_annotator'][confif_param] = param_value
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
        output_dir   = self._create_output_dirs()['result']
        inputs       = self._get_input_fpaths()
        repr_paths   = inputs['repr']


        output_files = {'vcf': ('vcf',),
                        'orf':('orf_seq.fasta', 'orf_pep.fasta'),
                        'ssr':('ssr',),
                        'gff':('gff',),
                        'orthologs':('orthologs',)}

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

            io_fhands = {'in_seq': open(seq_path.last_version),
                         'outputs': outputs}
            seq_pipeline_runner(pipeline=None,
                                configuration=None,
                                io_fhands=io_fhands)


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
            name   = filter_data['name']
            if 'step_name' in filter_data:
                step_name = filter_data['step_name']
            else:
                step_name = name
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

            filter_config['step_name'] = step_name
            configuration[name] =  filter_config
        pipeline = self._get_pipeline_from_step_name(configuration)

        return pipeline, configuration

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
        for name in configuration:
            step_name = configuration[name]['step_name']
            #we have to remove the step name from the configuration because this
            #dict will be used as **kargs
            del configuration[name]['step_name']
            step = translator[step_name]
            step['name_in_config'] = name
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

        #first we need some blasts
        project_dir = self._project_settings['General_settings']['project_path']
        blasts = {}
        for input_ in inputs['input']:
            input_fpath = input_.last_version
            for database in description_databases:
                db_kind = blast_settings[database]['kind']
                if db_kind == 'nucl':
                    blast_program = 'tblastx'
                else:
                    blast_program = 'blastx'

                blastdb = blast_settings[database]['path']
                blast   = backbone_blast_runner(query_fpath=input_fpath,
                                                project_dir=project_dir,
                                                blast_program=blast_program,
                                                blast_db=blastdb,
                                                dbtype=db_kind)
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
        repr_dir = output_dirs['repr_dir']
        result_dir = output_dirs['result']
        blast_settings = self._project_settings['blast']

        annot_settings = self._project_settings['Annotation']
        go_settings = annot_settings['go_annotation']
        go_database = go_settings['blast_database']

        if 'create_dat_file' in go_settings:
            create_dat  = go_settings['create_dat_file']
        else:
            create_dat  = None

        if 'java_memory' in  go_settings:
            java_memory = go_settings['java_memory']
        else:
            java_memory = None

        if 'prop_fpath' in go_settings:
            prop_fpath = go_settings['prop_fpath']
        else:
            prop_fpath = None

        #first we need some blasts
        project_dir = self._project_settings['General_settings']['project_path']
        chop_big_xml, num_items = True, 1000
        blasts = {}
        for input_ in inputs['input']:
            db_kind = blast_settings[go_database]['kind']
            if db_kind == 'nucl':
                blast_program = 'tblastx'
            else:
                blast_program = 'blastx'

            blastdb = blast_settings[go_database]['path']
            input_fpath = input_.last_version
            blast   = backbone_blast_runner(query_fpath=input_fpath,
                                            project_dir=project_dir,
                                            blast_program=blast_program,
                                            blast_db=blastdb,
                                            dbtype=db_kind)


            if chop_big_xml:
                #chopped_blast = open('/tmp/blast_itemized.xml', 'w')
                chopped_blast = NamedTemporaryFile(suffix='.xml')
                for blast_parts in xml_itemize(blast, 'Iteration', num_items):
                    chopped_blast.write(blast_parts)
                chopped_blast.flush()
                blast = chopped_blast
            #print "Itemized done"
            #raw_input()

            blasts[input_fpath] = blast

        # prepare pipeline
        pipeline      = [annotate_gos]
        configuration = {}

        for input_ in  inputs['input']:
            input_fpath = input_.last_version
            if create_dat:
                dat_fpath = os.path.join(result_dir, input_.basename + '.b2g.dat')
            else:
                dat_fpath = None
            annot_fpath = os.path.join(result_dir, input_.basename + '.b2g.annot')

            step_config = {'blast': blasts[input_fpath],
                           'dat_fpath': dat_fpath,
                           'annot_fpath': annot_fpath,
                           'java_memory':java_memory,
                           'prop_fpath':prop_fpath}

            configuration[input_.basename] = {'annotate_gos': step_config}

        result = self._run_annotation(pipeline=pipeline,
                                    configuration=configuration,
                                    inputs=inputs,
                                    output_dir=repr_dir)
        #remove the temporal files
        for input_ in  inputs['input']:
            blasts[input_fpath].close()
        return result


DEFINITIONS = {
    'filter_snvs':
        {'inputs':{
            'repr':
                {'directory': 'annotation_repr',
                 'file_kinds': 'sequence_files'},
            },
         'outputs':{'repr_dir':{'directory': 'annotation_repr'}},
         'analyzer': SnvFilterAnalyzer},
    'annotate_snv':
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
    'write_annotation':
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
    'annotate_description':
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

    'annotate_microsatellite':
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
    'annotate_orf':
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
    'annotate_go':
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
