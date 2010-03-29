'''
A blast runner hardly tied to backbone folder structure
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

from os.path import join, exists, splitext, basename
from os import makedirs, symlink
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES)
from franklin.utils.cmd_utils import call
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.seqio_utils import seqio
from tempfile import NamedTemporaryFile

def _get_basename(fpath):
    'It returns the base name without path and extension'
    return splitext(basename(fpath))[0]

def _create_temp_fasta_file(fpath):
    'It creates a fasta file format temfile'
    fasta_fhand = NamedTemporaryFile(suffix='.fasta', mode='a')
    seqio(in_seq_fhand=open(fpath),
          out_seq_fhand=fasta_fhand,
          out_format='fasta')
    return fasta_fhand

def backbone_blast_runner(query_fpath, project_dir, blast_program,
                          blast_db=None, blast_db_seq=None, dbtype='nucl'):
    '''It returns the blast if the results doesn't exist'''
    if blast_db is None and blast_db_seq is None:
        raise RuntimeError('It needs a blast database or seqfile')

    #create a logger
    logger = logging.getLogger("franklin")

    query_basename = _get_basename(query_fpath)
    blast_dir = join(project_dir, BACKBONE_DIRECTORIES['blast_dir'])

    if blast_db:
        result_dir = join(blast_dir, query_basename, _get_basename(blast_db))
    else:
        result_dir = join(blast_dir, query_basename,
                          _get_basename(blast_db_seq))
    if not exists(result_dir):
        makedirs(result_dir)
    result_fpath = join(result_dir,
                        '%s.%s.xml' % (BACKBONE_BASENAMES['blast_basename'],
                                       blast_program))
    if exists(result_fpath):
        logger.info('Using the stored blast result %s' % result_fpath)
        return open(result_fpath)

    #the input file should be fasta
    fasta_query_fhand = None
    fasta_db_fhand = None
    if guess_seq_file_format(open(query_fpath)) != 'fasta':
        fasta_query_fhand = _create_temp_fasta_file(query_fpath)
        query_fpath = fasta_query_fhand.name

    #we have to create a database in BACKBONE_DIRECTORIES['blast_databases']
    if blast_db_seq:
        #the name should be the basename of the blast_db_seq
        db_dir = join(project_dir, BACKBONE_DIRECTORIES['blast_databases'])
        if not exists(db_dir):
            makedirs(db_dir)
        db_seq_fpath = join(db_dir, _get_basename(blast_db_seq))
        if not exists(db_seq_fpath):
            #which is the name of the new databae?
            blast_db_seq_format = guess_seq_file_format(open(blast_db_seq))
            if blast_db_seq_format == 'fasta':
                symlink(blast_db_seq, db_seq_fpath)
            else:
                seqio(in_seq_fhand=open(blast_db_seq),
                      out_seq_fhand=open(db_seq_fpath, 'w'),
                      out_format='fasta')
            logger.info('Formatting the database %s' % db_seq_fpath)
            makeblastdb(db_seq_fpath, dbtype=dbtype)
        blast_db = db_seq_fpath

    logger.info('Running the blast %s' % result_fpath)
    blast_runner(query_fpath, blast_db, blast_program, result_fpath)

    if fasta_query_fhand:
        fasta_query_fhand.close()
    if fasta_db_fhand:
        fasta_db_fhand.close()

    return open(result_fpath)

def blast_runner(seq_fpath, blast_db, blast_type, result_fpath):
    'It runs a blast giving a file and a database path'
    cmd = [blast_type, '-db', blast_db, '-num_alignments', '25',
           '-num_descriptions', '25', '-evalue', '0.0001', '-outfmt', '5',
           '-query', seq_fpath, '-out', result_fpath]
    call(cmd, raise_on_error=True)

def makeblastdb(seq_fpath, dbtype, outputdb=None):
    'It creates the blast db database'
    cmd = ['makeblastdb', '-in', seq_fpath, '-dbtype', dbtype]
    if outputdb is not None:
        cmd.extend(['-out', outputdb])
    call(cmd, raise_on_error=True)
