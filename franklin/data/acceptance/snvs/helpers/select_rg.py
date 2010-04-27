import sys
rg = set()
for line in open(sys.argv[1]):
    line=line.strip()
    if line.startswith('@'):
        continue
    for item in line.split():
        if "RG:Z" in item:
            rg.add(item.split(':')[2])
print "\n".join(list(rg))
