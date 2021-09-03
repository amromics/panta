# Check wheather clusters change or not after add new samples

import csv
import sys
file1 = sys.argv[1]
file2 = sys.argv[2]

old_clusters = []
for row in csv.reader(open(file1, 'r')):
    this_cluster = set()
    if row[0] == 'Gene':
        continue
    for cell in row[8:]:
        genes = cell.split('\t')
        for gene in genes:
            if gene == '':
                continue
            this_cluster.add(gene)
    old_clusters.append([row[3], this_cluster])

new_clusters = []
for row in csv.reader(open(file2, 'r')):
    this_cluster = set()
    if row[0] == 'Gene':
        continue
    for cell in row[8:]:
        genes = cell.split('\t')
        for gene in genes:
            if gene == '':
                continue
            this_cluster.add(gene)
    new_clusters.append([row[3], this_cluster])

same = 0
difference = 0
total = 0
for i in old_clusters:
    subset = False
    for j in new_clusters:
        if i[1].issubset(j[1]):
            subset = True
    

    ls = [i[0], ]
    if subset == True:
        same += 1
        ls.append('Same')
    else:
        difference +=1
        ls.append('Diff')
    total += 1
    print('\t'.join(ls))

print(str(same) + '\t' + str(difference) + '\t'+ str(total))

