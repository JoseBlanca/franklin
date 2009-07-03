'''
Created on 2009 mai 21

@author: peio
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

from xml.sax import make_parser
from xml.sax.handler import ContentHandler
import StringIO
from biolib.biolib_utils import xml_itemize
from biolib.biolib_cmd_utils import  call

class InterProHandler(ContentHandler):
    '''  It parses the iprscan output'''

    #pylint: disable-msg=W0231
    def __init__(self):
        ''' It initiates the parser'''
        self.proteins  = []
        self._protein = None
        self.interpro = None
        self.match    = None
        self.is_match = False

    #pylint: disable-msg=C0103
    def startElement(self, name, attrs):
        '''It check in each line the key of the bracket and do what it needs '''
        if name == 'protein':
            self._protein = {}
            self._protein['id']       = attrs.get('id', None)
            self._protein['length']   = attrs.get('length', None)
            self._protein['interpro'] = []
        if name == 'interpro':
            self.interpro = {}
            self.interpro['id']   = attrs.get('id', None)
            self.interpro['name'] = attrs.get('name', None)
            self.interpro['type'] = attrs.get('type', None)
            self.interpro['clasifictions'] = []
            self.interpro['matches']       = []
        if name == 'classification':
            go_id = attrs.get('id', None)
            self.interpro['clasifictions'].append(go_id)
        if name == 'match':
            self.is_match = True
            self.match = {}
            self.match['id'] = attrs.get('id', None)
            self.match['name'] = attrs.get('name', None)
            self.match['dbname'] = attrs.get('dbname', None)
        if name == 'location' and  self.is_match:
            self.match['start']    = attrs.get('start', None)
            self.match['end']      = attrs.get('end', None)
            self.match['score']    = attrs.get('score', None)
            self.match['status']   = attrs.get('status', None)
            self.match['evidence'] = attrs.get('evidence', None)
    #pylint: disable-msg=C0103
    #pylint: disable-msg=W0613
    def endElement(self, name):
        ''' Signals the end of each element and performs actions'''
        if name == 'protein':
            self.proteins.append(self._protein)
        if name == 'interpro':
            self._protein['interpro'].append(self.interpro)
        if name == 'match':
            self.interpro['matches'].append(self.match)
            self.is_match = False
    def get_protein(self):
        '''It returns a protein each time'''
        protein      = self.proteins[0]
        self.proteins = self.proteins[1:]
        return protein

def xml_iprscan_parser_iter(fhand):
    '''It parses an iprscan xml result file.

    It's an iterator that yields the different protein sections
    '''
    parser = make_parser()
    handler = InterProHandler()
    parser.setContentHandler(handler)
    for protein_section in xml_itemize(fhand, 'protein'):
        fhand = StringIO.StringIO(protein_section)
        parser.parse(fhand)
        yield handler.get_protein()

def xml_iprscan_parser(fhand):
    '''It parses an iprscan xml result file.

    It returns a list with one dict for every protein section.
    '''
    parser = make_parser()
    handler = InterProHandler()
    parser.setContentHandler(handler)
    parser.parse(fhand)
    return handler.proteins

def iprscan_run(in_fasta_fpath, out_fpath):
    ''' It runs iprscan and returns the result in a variable. It assumes that
     the iprscan binary is correctly instaled'''
    iprscan_bin = 'iprscan'

    cmd = [iprscan_bin, '-cli', '-i', in_fasta_fpath, '-o', out_fpath,
           '-goterm', '-iprlookup','-format', 'xml' ]
    #pylint: disable-msg=W0612
    stdout, stderr, retcode = call(cmd)
    if retcode:
        raise RuntimeError(stderr)




