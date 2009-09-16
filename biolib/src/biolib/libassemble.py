'''
Created on 15/09/2009

@author: peio
'''
import os
from biolib.biolib_seqio_utils import (guess_seq_file_format,
                                       quess_seq_type, seqio)
from biolib.biolib_cmd_utils import call
def check_cfile(fhand):
    'It checks and prepares the configurations looking to the cfile'
    errors = []
    fhand.seek(0)
    conf = eval(fhand.read())
    # check_working_dir
    if 'work_dir' not in conf:
        errors.append('Working dir not configured')
    elif not os.path.exists(conf['work_dir']):
        os.mkdir(os.path.abspath(conf['work_dir']))

    # reference
    conf, errors = _check_seq_definitions(conf, errors, 'reference')
    # reads
    conf, errors = _check_seq_definitions(conf, errors, 'reads')

    # Is the cfile OK?
    if errors:
        error_msg = ["Please correct this errors in your config file:\n"]
        for error in errors:
            error_msg.append('\t.-%s\n' % error)
        raise RuntimeError("".join(error_msg))
    return conf

def _check_seq_definitions(conf, errors, kind):
    'It checks seq definition'
    if kind not in conf:
        errors.append('%s not in config file' % kind)
    else:
        for i, dict_def in enumerate(conf[kind]):
            conf[kind][i], errors = _check_seq_definition(dict_def, errors, kind)
    return conf, errors

def _check_seq_definition(dictionary, errors, kind):
    'It checks if the seq definitions dictionary is correct'
    #defaults aligners:
    aligners = {'short_seqs':{'binary':'nucmer' , 'parameters': '-l 20 -c 20'},
                'long_seqs' :{'binary':'nucmer' , 'parameters': ''}}

    #check required fields
    for field in ['name', 'seq_fpath']:
        if field not in dictionary:
            errors.append('%s : required %s field' % (kind, field))

    #qual file
    if 'qual_fpath' not in dictionary:
        dictionary['qual_fpath'] = None

    #check optional
    # check if we can guess file format
    if 'format' not in dictionary and not errors:
        try:
            format = guess_seq_file_format(open(dictionary['seq_fpath']))
            dictionary['format'] = format
        except ValueError:
            errors.append('File Format not defined and not recognizable')

    # check the seq type and use the apropiated alignments
    if kind == 'reads' and 'seq_type' not in dictionary and not errors:
        seq_type = quess_seq_type(open(dictionary['seq_fpath']), dictionary['format'], 40)
        dictionary['seq_type'] = seq_type
        if 'aliger' not in dictionary:
            dictionary['aligner'] = aligners[seq_type]['binary']
            dictionary['aligner_parameters'] = aligners[seq_type]['parameters']

    return dictionary, errors

def load_seqs_in_bank(configuration):
    '''This function loads the data in AMOS bank directory/format. It takes
     the input files from configuration'''
    seqs_to_load = get_files_to_load_in_bank(configuration)
    work_dir     = configuration['work_dir']
    os.chdir(work_dir)

    cmd = ['tarchive2amos', '-o', 'message_file' ]
    cmd.extend(seqs_to_load)
    stdout, stderr, retcode = call(cmd)

    if retcode:
        raise RuntimeError


def get_files_to_load_in_bank(configuration):
    '''This functions get files to load in bank from configuration.

    If it is needed it converts the files to fasta'''
    seqs_to_load = []
    work_dir     = configuration['work_dir']
    for seqs_file in configuration['reference'] + configuration['reads'] :
        seq_fpath  = seqs_file['seq_fpath']
        qual_fpath = seqs_file['qual_fpath']

        if 'format' not in seqs_file:
            format = guess_seq_file_format(open(seq_fpath))
        else:
            format     = seqs_file['format']
        if format == 'fasta':
            seqs_to_load.append(seq_fpath)
        else:
            in_seq_fhand = open(seq_fpath)

            if qual_fpath == None:
                in_qual_fhand = None
            else:
                in_qual_fhand = open(qual_fpath)

            basename = '.'.join(os.path.basename(seq_fpath).split('.')[:-1])
            out_seq_fpath = os.path.join(work_dir, basename + '.fasta')
            out_seq_fhand = open(out_seq_fpath, 'w')
            if format == 'fasta' and qual_fpath == None:
                out_qual_fhand = None
            else:
                out_qual_fpath = os.path.join(work_dir, basename+'.qual.fasta')
                out_qual_fhand = open(out_qual_fpath, 'w')

            seqio(in_seq_fhand, in_qual_fhand, format,
                  out_seq_fhand, out_qual_fhand, 'fasta')
            seqs_to_load.append(out_seq_fhand.name)

    return seqs_to_load




