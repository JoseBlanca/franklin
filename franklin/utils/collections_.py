'''
Created on 04/09/2009

@author: jose
'''

# Copyright 2009 Jose Blanca, Peio Ziarsolo, COMAV-Univ. Politecnica Valencia
# This file is part of franklin.
# project is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# franklin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR  PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with franklin. If not, see <http://www.gnu.org/licenses/>.

#import pysparse
import tempfile

class FileCachedList(object):
    '''A list cached in a file.

    The main aim is to store really big lists in files with an iterator
    interface.
    '''
    def __init__(self, type_):
        '''It creates a FileCachedList.

        It requires the type of objects that will be stored.
        '''
        self._type = type_
        #the file to store the data
        self._wfhand = tempfile.NamedTemporaryFile()
        #the file to read the data
        self._rfhand = None

    def append(self, item):
        'It writes one element at the file end'
        self._wfhand.write('%s\n' % str(item))

    def extend(self, items):
        'It adds a bunch of items from a list or from a FileCachedList'
        if 'items' in dir(items):
            items = items.items()
        for item in items:
            self.append(item)

    def items(self):
        'It yields all items'
        self._wfhand.flush()
        rfhand = open(self._wfhand.name)
        for line in rfhand:
            if len(line) == 0:
                raise StopIteration
            yield self._type(line.strip())

class RequiredPosition(object):
    '''This class checks if the given cromosome and positions are in the given
     file.

    This file must be ordered'''
    def __init__(self, fhand):
        '''It initializes the class '''
        self._fhand = fhand
        self._fhand.seek(0)
        self._cache = (None, None)

    def __getitem__(self, index):
        'The get item method. It look the location in the file'
        def bigger_than_index(cache_index, index):
            'Is bigger??'
            return (cache_index[0] > index[0] or
                   (cache_index[0] == index[0] and cache_index[1] > index[1]))

        index = index[0], int(index[1])
        cache_index = self._cache
        if cache_index == index:
            return True
        elif bigger_than_index(cache_index, index):
            return False
        for line in self._fhand:
            line = line.strip()
            if not line:
                continue
            cromosome, location = line.split()
            line_index  = cromosome, int(location)
            self._cache = line_index

            # if the given location is in the file
            if index == line_index:
                return True
            #if the asked_index > index
            elif bigger_than_index(line_index, index):
                break
        return False

def _ref_name(reference):
    'It returns the name of the reference'
    if 'name' in dir(reference):
        return reference.name
    else:
        return reference

def item_context_iter(items, window=None):
    '''Given an iter with Locatable items it returns an item, context iter,

    The items in the iter should have a location property (numerical) and a
    reference porperty. An item located based on its reference and its location.
    The iter should be sorted by references and locations.
    This generator will yield every item and its surrounding items (context).
    If window is not given all items from every reference will be included in
    the context.
    '''
    context = []
    if window is not None:
        width = window / 2.0
    else:
        width = None
    current_item = None #the item that we're yielding now
    items_buffer = []    #the items to be returned
    context_buffer = []
    right_edge_item = None
    last_item_in_iter = False
    while True:
        #if we have consumed the right item in the last iteration we ask for
        #other one
        if right_edge_item is None:
            try:
                right_edge_item = items.next()
            except StopIteration:
                #we have to empty the buffers before finishing
                last_item_in_iter = True
        #the case or the empty iter
        if not right_edge_item and not items_buffer:
            raise StopIteration
        #if the right item is in the same reference as the current one we add it
        #to the buffers
        if (right_edge_item is not None
            and (current_item is None or
                 _ref_name(right_edge_item.reference) ==
                                           _ref_name(current_item.reference))):
            #a buffer with the items to the right of the current one
            items_buffer.append(right_edge_item)
            #a buffer with the items to be added to the context
            context_buffer.append(right_edge_item)
            #we have consumed this right item
            right_edge_item = None
        #do we have to yield an item? we have to fill up the context buffer
        #first.
        #if there is no width but we're at the last item of the reference
        if (last_item_in_iter or
            #if we're at a new reference
            (_ref_name(items_buffer[-1].reference) !=
                                       _ref_name(items_buffer[0].reference)) or
            #if we're yet to close to the start we don't return anything
            (width is not None and items_buffer[-1].location > width)):
            current_item = items_buffer.pop(0)

        if not current_item:
            continue
        current_location = current_item.location
        current_reference = _ref_name(current_item.reference)
        #which items do we have to add to the context in the right side?
        while True:
            #if there are still items that might be added to the reference
            #and the item to be added has the same reference as the current
            if (context_buffer and
                   _ref_name(context_buffer[0].reference) == current_reference):
                if((width is None) or
                   (width is not None and
                        context_buffer[0].location - current_location < width)):
                    #and the item is close enough to the current item
                    context.append(context_buffer.pop(0))
                else:
                    break
            else:
                break
        #purge the items from the context that are not close to the current item
        while True:
            if not context:
                break
            if ((width is not None and
                 current_location - context[0].location  > width) or
                 current_reference != _ref_name(context[0].reference)):
                context.pop(0)
            else:
                break
        yield current_item, context
        #if we have consumed the buffers for this reference we have to go to the
        #next reference
        current_item = None
        if not items_buffer:
            current_item = None
            #if there's no more item and nothing in the buffers we're done
            if last_item_in_iter:
                raise StopIteration

def list_pairs_iter(items):
    'Given a list it yields all pairs between the items'
    length = len(items)
    for i_index in range(length):
        for j_index in range(i_index +1, length):
            yield(items[i_index], items[j_index])

def list_consecutive_pairs_iter(items):
    'Given a list it yields all consecutive pairs (1,2, 2,3, 3,4)'
    first_time = True
    first_item, second_item = None, None
    for item in items:
        if first_time:
            first_item = item
            first_time = False
            continue
        second_item = item
        yield first_item, second_item
        first_item, second_item = second_item, first_item


