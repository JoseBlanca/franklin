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

from os.path import join, exists, splitext, basename, split, abspath
from os import makedirs
from franklin.backbone.specifications import (BACKBONE_DIRECTORIES,
                                              BACKBONE_BASENAMES)
from franklin.utils.cmd_utils import call
from franklin.seq.readers import guess_seq_file_format
from franklin.utils.seqio_utils import seqio
from tempfile import NamedTemporaryFile


def _get_basename(fpath):
    'It returns the base name without path and extension'
    return splitext(basename(fpath))[0]

def backbone_blast_runner(query_fpath, project_dir, blast_program,
                          blast_db=None, blast_db_seq=None, dbtype='nucl'):
    '''It returns the blast if the results doesn't exist'''
    if blast_db is None and blast_db_seq is None:
        raise RuntimeError('It needs a blast database or seqfile')

    query_basename = _get_basename(query_fpath)
    blast_dir = join(project_dir, BACKBONE_DIRECTORIES['blast_dir'])

    result_dir = join(blast_dir, query_basename, _get_basename(blast_db))
    if not exists(result_dir):
        makedirs(result_dir)
    result_fpath = join(result_dir, BACKBONE_BASENAMES['blast_result_file'])
    if exists(result_fpath):
        return open(result_fpath)

    def _create_temp_fasta_file(fpath):
        'It creates a fasta file format temfile'
        fasta_fhand = NamedTemporaryFile(suffix='.fasta', mode='a')
        seqio(in_seq_fhand=open(fpath),
              out_seq_fhand=fasta_fhand,
              out_format='fasta')
        return fasta_fhand

    #the input file should be fasta
    if guess_seq_file_format(query_fpath) != 'fasta':
        fasta_query_fhand = _create_temp_fasta_file(query_fpath)
        query_fpath = fasta_query_fhand.name
    if blast_db_seq and guess_seq_file_format(blast_db_seq) != 'fasta':
        fasta_db_fhand = _create_temp_fasta_file(blast_db_seq)
        blast_db_seq = fasta_db_fhand.name

    if blast_db is None:
        blastdb = _guess_blastdb_path(blast_db_seq, project_dir)
        makeblastdb(blast_db_seq, dbtype=dbtype, outputdb=blastdb)

    blast_runner(query_fpath, blastdb, blast_program, result_fpath)

    fasta_query_fhand.close()
    fasta_db_fhand.close()

    return open(result_fpath)

def blast_runner(seq_fpath, blastdb, blast_type, result_fpath):
    'It runs a blast giving a file and a database path'
    cmd = [blast_type, '-db', blastdb, '-num_alignments', '20',
           '-num_descriptions', '20', -'evalue', '0.0001', '-outfmt', '5',
           '-query', seq_fpath, '-out', result_fpath]
    call(cmd, raise_on_error=True)

def _guess_blastdb_path(blastseq, project_dir):
    'It guess the blastdb path taking into account the seqpath'
    if 'annotations/input' in blastseq:
        basename_ = basename(blastseq)
        blast_db_dir = join(project_dir,
                            BACKBONE_DIRECTORIES['blast_databases'])
        if not exists(blast_db_dir):
            makedirs(blast_db_dir)
        return join(BACKBONE_DIRECTORIES['blast_databases'], basename_)
    else:
        return split(abspath(blastseq))[0]

def makeblastdb(seq_fpath, dbtype, outputdb=None):
    'It creates the blast db database'
    cmd = ['makeblastdb' '-in', seq_fpath, '-dbtype', dbtype]
    if outputdb is not None:
        cmd.extend(['-out', outputdb])
    call(cmd, raise_on_error=True)





