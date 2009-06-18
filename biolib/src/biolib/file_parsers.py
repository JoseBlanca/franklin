'''
Created on 2009 eka 18

@author: peio
'''

def library_parser(fhand):
    ''' It parses the file and yields a dictionary for each  library'''
    library = {}
    for line in fhand:
        line = line.strip()
        if line.startswith('library_definition'):
            #a new library starts
            if library:    #there was a previous library
                yield library
                library = {} #for the following library
            continue
        elif not line:
            continue
        elif line.startswith('format-version'):
            continue
        #if we're here line should be key: value1,value2, value3
        key, values = line.split(':', 1)
        key = key.strip()
        values = [value.strip() for value in values.split(',')]
        #some values should be scalars not lists
        if key == 'name':
            values = values[0]
        if key == 'organism':
            genus, specie = values[0].split(' ')
            library['genus']  = genus.strip()
            library['specie'] = specie.strip()
        elif key == 'type':
            cvname, cvtermname = values[0].split(':')
            library['cvname']     = cvname.strip()
            library['cvtermname'] = cvtermname.strip()
        else:
            library[key] = values
    else:
        #the last library
        yield library
