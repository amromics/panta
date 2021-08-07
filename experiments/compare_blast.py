import os
import sys
path = '/home/ntanh1999/amromics/amromics/pan-genome/tests/'
file1 = path + sys.argv[1]
file2 = path + sys.argv[2]


set1 = set()
for line in open(file1, 'r'):
    line = line.rstrip()
    cells = line.split('\t')
    set1.add(cells[0])
    set1.add(cells[1])

print(len(set1))


set2 = set()
for line in open(file2, 'r'):
    line = line.rstrip()
    cells = line.split('\t')
    set2.add(cells[0])
    set2.add(cells[1])

print(len(set2))


print(len(set1.intersection(set2)))
print(len(set1.difference(set2)))
print(len(set2.difference(set1)))