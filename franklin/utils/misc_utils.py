'''
Created on 2009 api 30

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

import tempfile, shutil
import os, re, math, subprocess
from UserDict import DictMixin
import franklin

DATA_DIR = os.path.join(os.path.split(franklin.__path__[0])[0], 'franklin',
                         'data')

def avail_phymem():
    """Return the amount of physical memory available, in bytes."""
    f = open('/proc/meminfo', 'r')
    free = None
    _flag = False
    for line in f:
        if line.startswith('MemFree:'):
            free = int(line.split()[1]) * 1024
            break
    f.close()
    return free

def cached_phymem():
    """Return the amount of cached memory on the system, in bytes.
    This reflects the "cached" column of free command line utility.
    """
    f = open('/proc/meminfo', 'r')
    for line in f:
        if line.startswith('Cached:'):
            f.close()
            return int(line.split()[1]) * 1024

def total_free_mem():
    "It returns total free memory. Taking into account cached memory"
    return avail_phymem() + cached_phymem()

def get_num_threads(threads, limit_by_memory=None):
    """It returns num of threads to use in parallel.

    You can pass to the funaction the  memory you want to use each thread.
    It calculates the number of treads
    In megabytes
    """
    phisical_threads = os.sysconf('SC_NPROCESSORS_ONLN')
    if not threads:
        return 1
    elif type(threads) == type(0):
        return threads
    elif limit_by_memory is not None:
        total_free_mem_ = total_free_mem()
        threads = total_free_mem_ / (limit_by_memory * 1024 * 1024)
        if threads > os.sysconf('SC_NPROCESSORS_ONLN'):
            return phisical_threads
        else:
            return threads
    else:
        return phisical_threads

def float_lists_are_equal(list1, list2):
    'Given two lists it checks that all floats are equal'
    for num1, num2 in zip(list1, list2):
        assert floats_are_equal(num1, num2)

def floats_are_equal(num1, num2):
    'Given two numbers it returns True if they are similar'
    if num1 == 0.0:
        if num2 == 0.0:
            return True
        else:
            return False
    log1 = math.log(float(num1))
    log2 = math.log(float(num2))
    return abs(log1 - log2) < 0.01

class NamedTemporaryDir(object):
    '''This class creates temporary directories '''
    def __init__(self):
        '''It initiates the class.'''
        self._name = tempfile.mkdtemp()

    def get_name(self):
        'Returns path to the dict'
        return self._name
    name = property(get_name)
    def close(self):
        '''It removes the temp dir'''
        if os.path.exists(self._name):
            shutil.rmtree(self._name)

    def __del__(self):
        '''It removes de temp dir when instance is removed and the garbaje
        colector decides it'''
        self.close()


def _remove_atributes_to_tag(tag):
    '''It removees atributes to a xml tag '''
    mod_tag = "".join(tag).split(' ')
    if len(mod_tag) > 1:
        return mod_tag[0] + '>'
    else:
        return mod_tag[0]

def _get_xml_header(fhand, tag):
    '''It takes the header of the xml file '''
    fhand.seek(0, 2)
    end_file = fhand.tell()

    fhand.seek(0, 0)
    header = []
    current_tag = []
    listed_tag = '<' + tag + '>'
    while True:
        if end_file <= fhand.tell():
            raise ValueError('End Of File. Tag Not found')

        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            header.append(letter)
        if letter == '>':
            mod_tag = _remove_atributes_to_tag(current_tag)
            if listed_tag == mod_tag:
                return  ''.join(header)
            else:
                header.extend(current_tag)
                current_tag = []
            in_tag = False

def _get_xml_tail(fhand, tag):
    '''It takes the tail of the xml file '''
    in_tag = False
    tail = []
    current_tag = []
    fhand.seek(-1, 2)
    listed_tag = list('</' + tag + '>')
    listed_tag.reverse()
    while True:
        if fhand.tell() == 0:
            raise ValueError('Start Of File. Tag Not found')
        letter = fhand.read(1)
        if letter == '>':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        else:
            tail.append(letter)

        if letter == '<':
            if current_tag == listed_tag:
                tail.reverse()
                return "".join(tail)
            else:
                tail.extend(current_tag)
                current_tag = []
            in_tag = False
        fhand.seek(-2, 1)

def xml_itemize(fhand, tag, num_items=1):
    '''It takes a xml file and it chunks it by the given key. It adds header if
    exists to each of the pieces. It is a generator'''
    fhand = get_fhand(fhand)
    fhand.seek(0, 2)
    end_file = fhand.tell()

    header = _get_xml_header(fhand, tag)
    tail = _get_xml_tail(fhand, tag)
    section = []
    current_tag = []
    listed_tag_s = '<' + tag + '>'
    listed_tag_e = '</' + tag + '>'
    in_tag, in_section = False, False

    fhand.seek(0, 0)

    items_in_buffer = 0
    buffer = []
    while True:
        if end_file <= fhand.tell():
            break
        letter = fhand.read(1)
        if letter == '<':
            in_tag = True
        if in_tag:
            current_tag.append(letter)
        if in_section:
            section.append(letter)
        if letter == '>':
            if listed_tag_s == _remove_atributes_to_tag(current_tag):
                in_section = True
                section.extend(current_tag)
            elif listed_tag_e == _remove_atributes_to_tag(current_tag):
                items_in_buffer += 1
                buffer.extend(section)
                if items_in_buffer >= num_items:
                    yield  header + "".join(buffer) + tail
                    items_in_buffer = 0
                    buffer = []
                section = []
                in_section = False
            in_tag = False
            current_tag = []
    #is there any reamining buffer
    if buffer:
        yield  header + "".join(buffer) + tail

class VersionedPath(object):
    'It represents a set of versioned files as one'
    def __init__(self, path):
        'It inits the object'
        self._path = os.path.abspath(path)
        if os.path.isdir(path):
            self.directory = path
            self.basename, version, self.extension = None, None, None
        else:
            self.directory, self.basename, version, self.extension = self._get_info(self._path)

    def _get_original_path(self):
        'it return original path'
        return self._path

    original_path = property(_get_original_path)

    @staticmethod
    def _get_info(fpath):
        'It returns all information about a given path'
        directory, fname = os.path.split(fpath)
        basename, ext = os.path.splitext(fname)
        extension = ext.replace('.', '', 1)
        match = re.search('(.*?)(\.\d+)?$', basename)
        basename, version = match.groups()
        version = int(version.lstrip('.')) if version else None
        return directory, basename, version, extension

    def __str__(self):
        'It returns the string representation'
        return os.path.join(self.directory,
                            self.basename + '.' + self.extension)
    def __repr__(self):
        'It returns the string representation'
        return '%s(%s)' % (self.__class__.__name__ , str(self))

    def _get_last_version(self):
        'It returns the path for the newest version of the file'
        return os.path.join(self.directory, self._get_last_fname())
    last_version = property(_get_last_version)

    def _get_next_version(self):
        'It returns the path for the next to the newest version of the file'
        directory, basename, version, extension = self._get_info(self.last_version)
        version = version + 1 if version is not None else 0
        return os.path.join(directory,
                            self._build_fname(basename, version, extension))
    next_version = property(_get_next_version)

    def list_fpaths(self):
        'It returns the fpaths found in the directory'
        return [os.path.join(self.directory, fname) for fname in
                                                     os.listdir(self.directory)]
    def list_fnames(self):
        'It returns the fnames found in the directory'
        return os.listdir(self.directory)

    def list_fpaths_versioned(self):
        'It returns the fpaths found in the directory (one for each versioned file)'
        fnames = os.listdir(self.directory)

        versioned_fnames = {}
        for fname in fnames:
            basename, version, extension = self._get_info(fname)[1:]
            if not (basename, extension) in versioned_fnames:
                versioned_fnames[basename, extension] = version
            else:
                newest_version = versioned_fnames[basename, extension]
                if newest_version is None or newest_version < version:
                    versioned_fnames[basename, extension] = version
        paths = []
        for (basename, extension), version in versioned_fnames.items():
            fname = self._build_fname(basename, version, extension)
            fpath = os.path.join(self.directory, fname)
            paths.append(fpath)
        return paths

    def list_paths_versioned(self):
        'It returns the fpaths found in the directory (one for each versioned file)'
        return [self.__class__(fpath) for fpath in self.list_fpaths_versioned()]

    def _build_fname(self, basename, version, extension):
        'It returns a valid fname'
        if version is None:
            return basename + '.' + extension
        else:
            return basename + '.' + str(version) + '.' + extension

    def _get_last_fname(self):
        'It takes the last version of a versioned file'
        regex = re.compile('(' + self.basename + ')\.(\d*).?(' + self.extension + ')$')
        last_file = None
        last_version = None
        #what happens if there is nothing in the directory
        for file_ in os.listdir(self.directory):
            match = regex.match(file_)
            if match:
                group = match.group(2)
                if group:
                    version = int(group)
                    if version > last_version:
                        last_version = int(version)
                        last_file = file_
                else:
                    if last_version is None:
                        last_file = file_
        if last_file is None:
            return self._build_fname(self.basename, None, self.extension)
        else:
            return last_file

    def __eq__(self, path):
        'The equality between two paths'
        dir1, basename1, ver1, ext1 = self._get_info(self._path)
        dir2, basename2, ver2, ext2 = self._get_info(str(path))
        return dir1 == dir2 and basename1 == basename2 and ext1 == ext2

    def __hash__(self):
        'It returns the hash'
        return hash(str(self))

#orderedDict taken from: http://code.activestate.com/recipes/576693/
class OrderedDict(dict, DictMixin):

    def __init__(self, *args, **kwds):
        if len(args) > 1:
            raise TypeError('expected at most 1 arguments, got %d' % len(args))
        try:
            self.__end
        except AttributeError:
            self.clear()
        self.update(*args, **kwds)

    def clear(self):
        self.__end = end = []
        end += [None, end, end]         # sentinel node for doubly linked list
        self.__map = {}                 # key --> [key, prev, next]
        dict.clear(self)

    def __setitem__(self, key, value):
        if key not in self:
            end = self.__end
            curr = end[1]
            curr[2] = end[1] = self.__map[key] = [key, curr, end]
        dict.__setitem__(self, key, value)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        key, prev, next = self.__map.pop(key)
        prev[2] = next
        next[1] = prev

    def __iter__(self):
        end = self.__end
        curr = end[2]
        while curr is not end:
            yield curr[0]
            curr = curr[2]

    def __reversed__(self):
        end = self.__end
        curr = end[1]
        while curr is not end:
            yield curr[0]
            curr = curr[1]

    def popitem(self, last=True):
        if not self:
            raise KeyError('dictionary is empty')
        if last:
            key = reversed(self).next()
        else:
            key = iter(self).next()
        value = self.pop(key)
        return key, value

    def __reduce__(self):
        items = [[k, self[k]] for k in self]
        tmp = self.__map, self.__end
        del self.__map, self.__end
        inst_dict = vars(self).copy()
        self.__map, self.__end = tmp
        if inst_dict:
            return (self.__class__, (items,), inst_dict)
        return self.__class__, (items,)

    def keys(self):
        return list(self)

    setdefault = DictMixin.setdefault
    update = DictMixin.update
    pop = DictMixin.pop
    values = DictMixin.values
    items = DictMixin.items
    iterkeys = DictMixin.iterkeys
    itervalues = DictMixin.itervalues
    iteritems = DictMixin.iteritems

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, self.items())

    def copy(self):
        return self.__class__(self)

    @classmethod
    def fromkeys(cls, iterable, value=None):
        d = cls()
        for key in iterable:
            d[key] = value
        return d

    def __eq__(self, other):
        if isinstance(other, OrderedDict):
            return len(self) == len(other) and self.items() == other.items()
        return dict.__eq__(self, other)

    def __ne__(self, other):
        return not self == other

class DisposableFile(file):
    'A file that removes the file when closed'
    def close(self):
        'This close removes the file when called'
        file.close(self)
        os.remove(self.name)

def get_fhand(file_, writable=False):
    'Given an fhand or and fpath it returns an fhand'
    if isinstance(file_, basestring):
        mode = 'w' if  writable else 'r'
        file_ = open(file_, mode)
    return file_

def _common_base(path1, path2):
    'it return the common and uncommon part of both strings'
    common = []
    path1 = path1.split(os.sep)
    path2 = path2.split(os.sep)


    for index, dir_ in enumerate(path1):
        try:
            if dir_ == path2[index]:
                common.append(dir_)
            else:
                break
        except IndexError:
            break

    uncommon1 = path1[len(common) :]
    uncommon2 = path2[len(common) :]
    return uncommon1, uncommon2

def _rel_path(path1, path2):
    'it return relative paths from path2 to path1'
    uncommon1, uncommon2 =  _common_base(path1, path2)

    path =  ['..'] * (len(uncommon2) - 1)
    path.extend(uncommon1)
    return os.sep.join(path)

def _check_call(cmd):
    'It checks the call and if retcode raises a nice error'
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    stderr = process.communicate()[1]
    retcode = process.poll()
    if retcode:
        msg = 'Problem creating relative path, command was: %s\nerror was: %s\n'
        msg %= ' '.join(cmd), stderr
        raise ValueError(msg)
    return

def rel_symlink(path1, path2):
    'It makes the relative symlink'
    path1 = os.path.abspath(path1)
    path2 = os.path.abspath(path2)
    dest_dir, fname2 = os.path.split(path2)
    rel_path1 = _rel_path(path1, path2)

    #we need a temp dir to be thread safe when creating the link
    #we cannot use chdir, that is not threadsafe
    temp_dir = tempfile.mkdtemp()
    temp_link_path = os.path.join(temp_dir, fname2)

    cmd = ['ln', '-s', '-T', rel_path1, temp_link_path]
    _check_call(cmd)
    cmd2 = ['mv', temp_link_path, path2]
    _check_call(cmd2)

    os.rmdir(temp_dir)

def get_db_connection(db_data):
    'It connects to a db using python dbi api. Returns a connection'
    if db_data['pass'] is None:
        db_data['pass'] = ''
    if db_data['adaptor'] == 'mysql':
        import MySQLdb
        conn = MySQLdb.connect(db=db_data['name'],
                               user=db_data['user'],
                               host=db_data['host'],
                               passwd=db_data['pass'])
    elif db_data['adaptor'] == 'psql':
        import psycopg2
        conn = psycopg2.connect(database = db_data['name'],
                                user = db_data['user'],
                                password = db_data['pass'],
                                host = db_data['host'])


    return conn


