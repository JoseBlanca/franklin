import sys

unigenes = [unigene.strip() for unigene in open('used_unigenes')]

for line in open(sys.argv[1]):
    if line.split()[2] in unigenes:
        print line.strip()



