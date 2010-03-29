'''
Created on 15/03/2010

@author: peio
'''
from tempfile import NamedTemporaryFile
import os, shutil
from franklin.utils.misc_utils import NamedTemporaryDir
from franklin.backbone.analysis import (Analyzer, scrape_info_from_fname,
                                        LastAnalysisAnalyzer)
from franklin.utils.cmd_utils import call
from franklin.seq.readers import guess_seq_file_format
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES)
from franklin.utils.seqio_utils import (seqio, cat, seqs_in_file,
                                        write_seqs_in_file)

class PrepareWSGAssemblyAnalyzer(Analyzer):
    '''It collects the cleaned reads to use by wsg. Wsg only uses reads with
    quality, so be will dismiss these sequences'''

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        self._log({'analysis_started':True})
        # TODO
        # Now we need to convert all the fastq in fasta and qual. If fasta
        #without qual, we need to give it a quality from settings. It can change
        # if people in celera get fastqToCA to work.
        tempdir = NamedTemporaryDir()
        tempdir_name = tempdir.name
        fasta_qual_files = []
        for seqfile in self._get_input_fpaths()['reads']:
            seqfhand = open(seqfile)
            basename = os.path.splitext(os.path.basename(seqfile))[0]
            fasta_fhand = open(os.path.join(tempdir_name,
                                            basename + '.fasta'), 'w')
            qual_fhand  = open(os.path.join(tempdir_name,
                                            basename + '.qual'), 'w')
            fasta_qual_files.append((fasta_fhand.name, qual_fhand.name))

            file_format = guess_seq_file_format(seqfhand)
            seqs = seqs_in_file(seqfhand, format=file_format)

            write_seqs_in_file(seqs,  seq_fhand=fasta_fhand,
                               qual_fhand=qual_fhand,
                               format='fasta')
            fasta_fhand.close()
            qual_fhand.close()
            seqfhand.close()

        # once we have the fasta and qual files, we need to convert them to WSG
        # format FRG. For fasta convertTo

        for fasta_fpath, qual_fpath in fasta_qual_files:
            file_info = scrape_info_from_fname(fasta_fpath)
            library   = file_info['lb']
            platform  = file_info['pl']

            cmd = ['convert-fasta-to-v2.pl', '-noobt', '-l', library, '-s',
                   fasta_fpath, '-q', qual_fpath]

            if platform == '454':
                cmd.append('-454')
            stdout = call(cmd, raise_on_error=True)[0]
            basename  = os.path.splitext(os.path.basename(fasta_fpath))[0]
            frg_fpath = os.path.join(tempdir.name, basename + '.frg')
            frg_fhand = open(frg_fpath, 'w')
            frg_fhand.write(stdout)
            frg_fhand.close()

        # Finally we do a cat of all the frg files ant write it to real output
        # dir
        output_dir = self._create_output_dirs()['assembly_input']
        frg_fhands = []
        for frg in os.listdir(tempdir.name):
            if frg.endswith('.frg'):
                frg_fhands.append(open(os.path.join(tempdir.name, frg)))
        final_frg = os.path.join(output_dir, BACKBONE_BASENAMES['merged_frg'])
        cat(frg_fhands, open(final_frg, 'w'))

        #close all files
        for frg_fhand in frg_fhands:
            frg_fhand.close()
        tempdir.close()
        self._log({'analysis_finished':True})

class WSGAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads using WSG assembler'

    def run(self):
        '''It runs the analysis.'''
        self._log({'analysis_started':True})
        proj_name = self._get_project_name()
        #we need an assembly dir for this run
        #in this case every assembly is done in a new directory
        #the directory for this assembly
        output_dirs = self._create_output_dirs(timestamped=True)
        assembly_dir = output_dirs['analysis']
        wsg_dir  = os.path.join(assembly_dir, '%s_assembly' % proj_name)

        frgs = self._get_input_fpaths()['frgs']

        stdout = open(os.path.join(assembly_dir, 'stdout.txt'), 'w')
        stderr = open(os.path.join(assembly_dir, 'stderr.txt'), 'w')
        cmd = ['runCA', '-d', wsg_dir, '-p', proj_name]
        cmd.extend(frgs)
        #print ' '.join(cmd)
        raw_input()
        retcode = call(cmd, stdout=stdout, stderr=stderr)[-1]
        if retcode:
            stdout.flush()
            stdout.seek(0)
            raise RuntimeError(stdout.read())
        stdout.close()
        stderr.close()

        # symlinks result to result dir
        fname = '%s.asm' % proj_name
        asm_result_fpath = os.path.join(output_dirs['analysis'], fname)
        os.symlink(asm_result_fpath, os.path.join(output_dirs['result'], fname))
        self._log({'analysis_finished':True})

class PrepareMiraAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads'

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        self._log({'analysis_started':True})
        files_454 = []
        files_sanger_with_qual = []
        files_sanger_without_qual = []
        for fpath in self._get_input_fpaths()['reads']:
            fhand = open(fpath)
            fname = os.path.split(fpath)[-1]
            if 'pl_454' in fname.lower():
                files_454.append(fhand)
            elif 'pl_sanger' in fname.lower():
                format_ = guess_seq_file_format(fhand)
                if format_ == 'fasta':
                    files_sanger_without_qual.append(fhand)
                elif format_ == 'fastq':
                    files_sanger_with_qual.append(fhand)

        #fastq are processed before
        files_sanger = files_sanger_with_qual[:]
        files_sanger.extend(files_sanger_without_qual)

        #all files should be fasta and fasta.qual
        output_dir = self._create_output_dirs()['assembly_input']
        project_name = self._get_project_name()
        for ext, files in (('_in.454', files_454),
                           ('_in.sanger', files_sanger)):
            base_name = os.path.join(output_dir, project_name + ext)
            fasta_fpath = base_name + '.fasta'
            qual_fpath = base_name + '.fasta.qual'
            if os.path.exists(fasta_fpath) or not files:
                continue
            fasta_fhand = open(fasta_fpath, 'w')
            qual_fhand = open(qual_fpath, 'w')
            self._cat_to_fasta(files, fasta_fhand, qual_fhand)
            fasta_fhand.close()
            qual_fhand.close()

        # close all files
        for file_ in files_454 + files_sanger:
            file_.close()
        self._log({'analysis_finished':True})

    @staticmethod
    def _files_to_temp_fasta(files):
        'It converts the given files to a temporary fasta and qual'
        fastas, quals = [], []
        for file_ in files:
            #are we dealing with a fastq file (with qual)
            if 'fastq' in os.path.splitext(file_.name)[-1]:
                qual  = NamedTemporaryFile(suffix='.qual')
                fasta = NamedTemporaryFile(suffix='.fasta')
                seqio(in_seq_fhand=file_, out_seq_fhand=fasta,
                        out_qual_fhand=qual, out_format='fasta')
            else:
                #the file is already fasta
                fasta = file_
                qual = None
            fastas.append(fasta)
            quals.append(qual)
        return fastas, quals

    def _cat_to_fasta(self, files, fasta_fhand, qual_fhand):
        'It cats the given files together and it returns a fasta and qual'
        #all files should be in fasta and qual
        fastas, quals = self._files_to_temp_fasta(files)
        cat(fastas, fasta_fhand)
        cat(quals , qual_fhand)

class MiraAssemblyAnalyzer(Analyzer):
    'It assembles the cleaned reads'

    def run(self):
        '''It runs the analysis. It checks if the analysis is already done per
        input file'''
        self._log({'analysis_started':True})
        proj_name = self._get_project_name()
        #we need an assembly dir for this run
        #in this case every assembly is done in a new directory
        #the directory for this assembly
        output_dirs = self._create_output_dirs(timestamped=True)
        assembly_dir = output_dirs['analysis']
        mira_dir  = os.path.join(assembly_dir, '%s_assembly' % proj_name)

        os.chdir(assembly_dir)
        #with technologies have we used?
        #each input file should have a link in the assembly dir
        techs = set()
        for fpath in self._get_input_fpaths()['reads']:
            fname = os.path.split(fpath)[-1]
            mira_input_fname = os.path.join(assembly_dir, fname)
            os.symlink(fpath, mira_input_fname)
            if '454' in fname:
                techs.add('454')
            elif 'sanger' in fname:
                techs.add('sanger')

        settings = self._project_settings['Mira']
        #job part of the command

        job = settings['job_options']
        job.extend(techs)
        job = ','.join(job)
        job = '-job=' + job

        cmd = ['mira', '-project=%s' % proj_name, '-fasta', job]
        if 'general_settings' in settings:
            general_settings = settings['general_settings']
            cmd.extend(general_settings)

        cmd.append('-OUT:rld=yes')

        for tech in techs:
            tech_str = '%s_settings' % tech
            if tech_str in settings:
                cmd.append(tech_str.upper())
                cmd.extend(settings[tech_str])

        stdout = open(os.path.join(assembly_dir, 'stdout.txt'), 'w')
        stderr = open(os.path.join(assembly_dir, 'stderr.txt'), 'w')
        retcode = call(cmd, stdout=stdout, stderr=stderr)[-1]
        if retcode:
            stdout.flush()
            stdout.seek(0)
            raise RuntimeError(stdout.read())

        stdout.close()
        stderr.close()
        #remove the log directory
        mirachkpt_dir = os.path.join(mira_dir,  "%s_d_chkpt" % proj_name)
        if os.path.exists(mirachkpt_dir):
            shutil.rmtree(mirachkpt_dir)

        #results
        mira_results_dir = os.path.join(mira_dir, '%s_d_results' % proj_name)
        results_dir = output_dirs['result']
        mira_basename = '%s_out' % proj_name
        file_exts = (('.ace', '.ace'),
                     ('.caf', '.caf'),
                     ('.unpadded.fasta', '.fasta'),
                     ('.unpadded.fasta.qual', '.qual'))
        for mira_ext, result_ext in file_exts:
            mira_fpath = os.path.join(mira_results_dir,
                                      mira_basename + mira_ext)
            res_fpath = os.path.join(results_dir,
                                     BACKBONE_BASENAMES['contigs'] + result_ext)
            if os.path.exists(mira_fpath):
                os.symlink(mira_fpath, res_fpath)


        mira_info_dir = os.path.join(mira_dir, '%s_d_info' % proj_name)
        if os.path.exists(mira_info_dir):
            os.symlink(mira_info_dir, os.path.join(results_dir,
                                                  BACKBONE_DIRECTORIES['info']))
        self._log({'analysis_finished':True})

DEFINITIONS ={
    'prepare_wsg_assembly':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'assembly_input':{'directory': 'assembly_input'}},
         'analyzer': PrepareWSGAssemblyAnalyzer,
        },
    'wsg_assembly':
        {'inputs':{
            'frgs':
                {'directory': 'assembly_input',
                 'file_kinds': 'frg_file'}
            },
         'outputs':{
                    'analysis':  {'directory': 'assemblies'},
                    'result':  {'directory': 'assembly_result'},
                    },
         'analyzer': WSGAssemblyAnalyzer,
        },
    'prepare_mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'cleaned_reads',
                 'file_kinds': 'sequence_files'}
            },
         'outputs':{'assembly_input':{'directory': 'assembly_input'}},
         'analyzer': PrepareMiraAssemblyAnalyzer,
        },
    'mira_assembly':
        {'inputs':{
            'reads':
                {'directory': 'assembly_input',
                 'file_kinds': 'sequence_or_qual_files'}
            },
         'outputs':{
                    'analysis':  {'directory': 'assemblies'},
                    'result':  {'directory': 'assembly_result'},
                    },
         'analyzer': MiraAssemblyAnalyzer,
        },
    'select_last_assembly':
        {'inputs':{'analyses_dir':{'directory': 'assemblies'}},
         'outputs':{'result':{'directory': 'assembly_result',
                              'create': False}},
         'analyzer': LastAnalysisAnalyzer,
        },
   }
