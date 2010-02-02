'It gives unique names for the features in the databases'

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


import sqlalchemy, re, os
from sqlalchemy import (Table, Column, Integer, String, Boolean, ForeignKey,
                        DateTime)
from biolib.db.db_utils import setup_mapping

from biolib.utils.seqio_utils import seqs_in_file, write_seqs_in_file
from datetime import datetime

def create_naming_database(engine):
    'It creates a new empty database to hold the naming schema status'
    #the table definition
    metadata = sqlalchemy.MetaData()
    metadata.bind = engine
    Table('names', metadata,
         Column('id', Integer, primary_key=True),
         Column('project_id', String(50), ForeignKey("projects.id"),
                nullable=False),
         Column('name', String(50), unique=True, nullable=False),
         Column('feature_type', String(50), nullable=False),
         Column('lock', Boolean),
         Column('description', String(200)),
         Column('date', DateTime)
    )
    Table('projects', metadata,
         Column('id', Integer, primary_key=True),
         Column('short_name', String(50), unique=True, nullable=False),
         Column('code', String(2), unique=True, nullable=False),
         Column('description', String(300))
    )
    metadata.create_all(engine)

def _setup_naming_database_mapping(engine):
    'It creates the orm mapping form the naming db'
    mapping_definitions = [
                    {'name':'projects'},
                    {'name':'names',
                     'relations':{'project_id':{'kind':'one2many',
                                                'rel_attr':'project'}}},
                    ]
    return setup_mapping(engine, mapping_definitions)


def add_project_to_naming_database(engine, name, code, description=None):
    'It adds a new project to the naming database'
    #the mapping
    #pylint: disable-msg=W0612
    row_classes = _setup_naming_database_mapping(engine)[1]
    session_klass = sqlalchemy.orm.sessionmaker(bind=engine)
    session = session_klass()
    project = row_classes['projects']()
    project.short_name = name
    project.code       = code
    if description is not None:
        project.description = description
    session.add(project)
    session.commit()

def project_in_database(engine, name):
    'It checks if the project is already added to the database'
    row_classes = _setup_naming_database_mapping(engine)[1]
    session_klass = sqlalchemy.orm.sessionmaker(bind=engine)
    session = session_klass()
    project_klass = row_classes['projects']
    try:
        session.query(project_klass).filter_by(short_name=name).one()
    except sqlalchemy.orm.exc.NoResultFound:
        return False
    return True

class _CodeGenerator(object):
    '''This class gives the next code, giving the last one'''
    # pylint: disable-msg=R0903
    #no need for more public methods
    def __init__(self, seed):
        '''init'''
        self._last_code  = seed
        self._dictionary = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    def next(self):
        '''It return the next code '''
        last_code = self._last_code
        length    = len(last_code)
        letter_part, number_part = '', ''
        next_code = None
        for item  in last_code:
            if item.isdigit():
                number_part += item
            else:
                letter_part += item
        len_letter, len_number = len(letter_part), len(number_part)
        if number_part and number_part == len_number * '9':
            if len_letter == 0:
                next_code = "A".ljust(length, "0")
            else:
                next_code = self._next_letter_code(letter_part).ljust(length,
                                                                          "0")
        elif len_number == 0:
            if letter_part == length * "Z":
                next_code = "0" * (length + 1)
            else:
                next_code = self._next_letter_code(letter_part)
        else:
            letters = str(int(number_part) + 1).rjust(len_number, '0')
            next_code = letter_part + letters
        self._last_code = next_code
        return next_code
    def _next_letter(self, letter):
        '''It returns the next letter one caracter'''
        pos = self._dictionary.index(letter)
        return self._dictionary[pos + 1]
    def _next_letter_code(self, letter):
        '''It return the nextletter code'''
        length = len(letter)
        if letter ==  length * "Z":
            return "A" * (length + 1)
        elif letter[-1] == "Z":
            return self._next_letter_code(letter[:-1]) + "A"
        else:
            return letter[:-1] + self._next_letter(letter[-1])

class DbNamingSchema(object):
    'This class gives unique names for the objects in the databases'
    def __init__(self, engine, project=None, feature_kind=None):
        '''It inits the NamingSchema.

        engine is a sqlalchemy engine db.
        project is the name of the project. Independent names will be created
        for each project.
        '''
        #the session
        session_klass = sqlalchemy.orm.sessionmaker(bind=engine)
        self._session = session_klass()

        self._code_generators = {}
        self._curr_code = {}  #the last code given for every feature_kind
        #pylint: disable-msg=W0612
        self._row_classes = _setup_naming_database_mapping(engine)[1]
        #the project instance
        self.project = project
        self.kind = feature_kind

    def set_project(self, name):
        'It sets the project property'
        if name is None:
            self._project = None
        else:
            project_klass = self._row_classes['projects']
            self._project = \
             self._session.query(project_klass).filter_by(short_name=name).one()
    def get_project(self):
        'It returns the project name'
        if self._project is None:
            return None
        else:
            return self._project.short_name
    project = property(get_project, set_project)

    def _get_so_term(self):
        'It returns the Sequence Ontology term  for the given feature'
        so_terms = {'EST':'SO:0000345', 'transcribed_cluster':'SO:0001457',
                    'sa':'SO:0001059'
             }
        if not self.kind in so_terms:
            raise ValueError('Unkown feature_kind, please use so')
        return so_terms[self.kind]

    def _get_type_code(self):
        'It returns the two letter code for the feature kind.'
        type_code = {'SO:0000345': 'ES', 'transcribed_cluster':'TC',
                     'EST':'ES', 'library':'lc', 'sa':'sa'}
        if self.kind not in type_code:
            raise ValueError('No code to the feature kind: ' + self.kind)
        return type_code[self.kind]

    def _get_code_generator(self):
        'It returns a code generator for the given object kind'
        if self.kind in self._code_generators:
            return self._code_generators[self.kind]
        #which was the last name given for this kind of feature
        names_klass = self._row_classes['names']
        last_name = \
            self._session.query(names_klass).filter_by(project=self._project,
               feature_type=self.kind).order_by(names_klass.date.desc()).first()
        if last_name:
            last_name = last_name.name
        else:
            #it's the first time we're using this type in this database
            #we add one row to the last_name table
            last_name = self._project.code + self._get_type_code() + '000000'
        code_gen = _CodeGenerator(last_name[4:])
        self._code_generators[self.kind] = code_gen
        return code_gen

    def get_uniquename(self, name=None):
        'It returns the valid name for the next feature in the project'
        kind = self.kind
        #pylint: disable-msg=W0613
        code_gen = self._get_code_generator()

        #now we want the next number
        num_next = code_gen.next()
        next_code = self._project.code + self._get_type_code() + num_next
        self._curr_code[kind] = next_code
        return next_code

    def commit(self, description=None):
        '''It stores the current code in the database.

        The next time that someone ask for a code to the database, this one
        will be used as seed.
        '''
        names_klass = self._row_classes['names']
        for kind, last_name in self._curr_code.items():
            new_row = names_klass(project=self._project,
                                  name=last_name,
                                  feature_type=kind,
                                  date=datetime.now(),
                                  description=description)
            self._session.add(new_row)
        self._session.commit()

    def get_names_from_db(self):
        'It shows all the names in the database'
        names_klass = self._row_classes['names']
        query = self._session.query(names_klass)
        if self.project is not None:
            query = query.filter_by(project=self._project)
        if self.kind is not None:
            query = query.filter_by(feature_type=self.kind)

        for name in query.all():
            name_dict = {}
            for field in ('name', 'feature_type', 'description', 'date'):
                name_dict[field] = getattr(name, field)
            name_dict['project'] = getattr(name, 'project').short_name
            yield name_dict

    def revert_last_name(self):
        '''It reverts the database to the previous state for a certain kind and
        project'''
        if self.project is None or self.kind is None:
            msg = 'Project and kind should be given to remove a name'
            raise ValueError(msg)
        names_klass   = self._row_classes['names']
        last_name = \
            self._session.query(names_klass).filter_by(project=self._project,
               feature_type=self.kind).order_by(names_klass.date.desc()).first()
        self._session.delete(last_name)
        self._session.commit()

NAMES_RE = {'fasta'  :{'sequence_names':[r'^>([^ \n]+).*$']},
            'library':{'library_names' :[r'\s*name\s*:\s*(\w+)']},
            'ace'    :{'contig_names':[],
                       'read_names'  :[]},
            'caf'    :{'contig_names':[],
                       'read_names'  :[]}
          }
class FileNamingSchema(object):
    '''It takes a naming file and it converts to a dict '''
    def __init__(self, fhand, naming_schema=None, feature_kind=None):
        '''The initiator '''
        self._fhand           = fhand
        self._naming_schema   = naming_schema
        self._naming_dict     = {}
        self._new_naming_dict = {}
        self._naming_file_to_dict()
        self._feature_kind = feature_kind

    def _set_kind(self, kind):
        'It sets the kind in of feature'
        if self._naming_schema:
            self._naming_schema.kind = kind
    def _get_kind(self):
        'It returns the kind of feature'
        if self._naming_schema:
            return self._naming_schema.kind
        else:
            return self._feature_kind
    kind = property(_get_kind, _set_kind)

    def get_uniquename(self, name=None):
        '''Given a name and a it returns a uniquename'''
        try:
            uniquename = self._naming_dict[name]
        except KeyError:
            try:
                uniquename = self._new_naming_dict[name]
            except KeyError:
                uniquename = None
        if uniquename is None and self.kind is None:
            msg = 'You must provide feature type if the name is not added'
            raise ValueError(msg)
        if uniquename is None:
            if self._naming_schema is None:
                msg = 'Uncached name and no naming schema'
                raise ValueError(msg)
            else:
                uniquename = self._naming_schema.get_uniquename()
                self._new_naming_dict[name] = uniquename
        return uniquename

    def commit(self):
        '''It commits the new  '''
        self._naming_schema.commit()
        self._write_dict_to_file()

    def _write_dict_to_file(self):
        '''It writes the new uniquenames to the fhand '''
        if not os.path.exists(self._fhand.name):
            mode = 'w'
        else:
            mode = 'a'
        self._fhand = open(self._fhand.name, mode)
        for name, uniquename in self._new_naming_dict.items():
            self._fhand.write('%s:%s\n' % (name, uniquename))
            self._new_naming_dict = {}
        self._fhand.flush()

    def _naming_file_to_dict(self):
        '''Giving a a file with name: uniquename translation it returns a
        dictionary'''
        for line in open(self._fhand.name, 'r'):
            if not line.isspace():
                items      = line.split(':')
                name       = items[0].strip()
                uniquename = items[1].strip()
                self._naming_dict[name] = uniquename

REPLACE_RE = {
    'fasta':[(r'^(>)([^ \n]+)(.*)$', 2)],
    'caf'  :[(r'^(DNA *: *)([^ \n]+)(.*)$', 2),
             (r'^(BaseQuality *: *)([^ \n]+)(.*)$', 2),
             (r'^(Sequence *: *)([^ \n]+)(.*)$', 2),
             (r'^(Assembled_from *)([^ ]+)( .*)$', 2)],
    'ace'  :[(r'^(CO +)([^ ]+)(.*)$', 2),
             (r'^(AF +)([^ ]+)(.*)$', 2),
             (r'^(RD +)([^ ]+)(.*)$', 2)]
}

def _change_names_in_files_regex(fhand_in, fhand_out, naming, file_format):
    '''It changes the accession names in the files.

    Given a naming schema and a file kind
    names found in the in file and it writes the result in the output file.
    '''
    def repl_factory(regex, group_num):
        '''It creates a repl function for the re.sub function.

        We need one of these function for each regular expression given.
        Besides the regex we need the group number that has the name to be
        changed.
        '''
        #the function
        def repl_fun(matchobj):
            'It returns the replacement for the match.group(1) string'
            replace = []
            for group_i, group in enumerate(matchobj.groups()):
                if group_i + 1 == group_num:
                    replace.append(naming.get_uniquename(name=group))
                else:
                    replace.append(group)
            return ''.join(replace)
        #the re object
        re_obj = re.compile(regex)
        return repl_fun, re_obj

    #the repl and re_obj for each regular expression
    re_list = REPLACE_RE[file_format]
    regex_list = []
    for regex in re_list:
        func, re_obj = repl_factory(regex[0], regex[1])
        regex_list.append({'function':func, 're':re_obj})

    fhand_in.seek(0)
    fhand_out.seek(0)
    for line in fhand_in:
        for regex in regex_list:
            line = re.sub(regex['re'], regex['function'], line)
        fhand_out.write(line)

def _change_names_in_files_by_seq(fhand_in, fhand_out, naming, file_format):
    'It replaces the seq name using the  per_seq method'
    seqs = seqs_in_file(fhand_in, format=file_format)

    for seq in seqs:
        new_name = naming.get_uniquename()
        seq.name = new_name
        seq.id   = new_name
        write_seqs_in_file([seq], fhand_out, format=file_format)

def change_names_in_files(fhand_in, fhand_out, naming, file_format):
    '''It changes the accession names in the files.

    Given a naming schema and a file kind
    names found in the in file and it writes the result in the output file.
    '''
    if file_format in REPLACE_RE:
        _change_names_in_files_regex(fhand_in, fhand_out, naming, file_format)
    else:
        _change_names_in_files_by_seq(fhand_in, fhand_out, naming, file_format)

#pylint: disable-msg=R0903
class _GeneralNameParser(object):
    '''This class parses some kind of files looking for names '''
    def __init__(self, fhand, kind):
        '''Initiator'''
        self.fhand  = fhand
        self.kind   = kind
        self.regexs = NAMES_RE[self.kind]
        self._names = {}
        self._parser()
        self._setup_accesor_methods()

    def _parser(self):
        '''The real parser '''
        for obj_kind in self.regexs:
            self._names[obj_kind] = []
        for line in self.fhand:
            for obj_kind, regex_list in self.regexs.items():
                for regex in regex_list:
                    match_obj = re.match(regex, line)
                    try:
                        name      = match_obj.groups()[0]
                        self._names[obj_kind].append(name)
                    #pylint: disable-msg=W0704
                    except AttributeError:
                        pass
    def _setup_accesor_methods(self):
        '''It creates classes for all kinds '''
        for obj_kind in self.regexs:
            def accessor():
                'It creates a method for each kind'
                return self._names[obj_kind]
            setattr(self, obj_kind, accessor)
