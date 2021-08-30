# so sánh pipeline vs Roary

import csv
import sys
path = '/home/ntanh1999/amromics/amromics/pan-genome/tests/'
file1 = path + sys.argv[1]+'/gene_presence_absence.csv'
file2 = path + sys.argv[2]+'/gene_presence_absence.csv'

list1 = []
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
    list1.append(this_cluster)

list2 = []
for row in csv.reader(open(file2, 'r')):
    this_cluster = set()
    if row[0] == 'Gene':
        continue
    for cell in row[14:]:
        genes = cell.split('\t')
        for gene in genes:
            if gene == '':
                continue
            this_cluster.add(gene)
    list2.append(this_cluster)

giong = 0
khac = 0
for i in list2:
    maxx = 0
    match = None
    
    for j in list1:
        intersection = i.intersection(j)
        if len(intersection) > maxx:
            maxx = len(intersection)
            match = j

    if match == None:
        khac += 1
        continue    
    chung = i.intersection(match)
    khac_1 = i.difference(match)
    khac_2 = match.difference(i)

    if len(khac_1) == 0 and len(khac_2) == 0:
        giong += 1
    else:
        khac += 1

tong = giong + khac
print(str(giong) + '\t' + str(khac) + '\t'+ str(tong))


# python3 compare_result.py Sa110_add_new_10 Sa110_roary_nosplit