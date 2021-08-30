# compare clusters in more details 

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

chap_nhan = 0
khong_chap_nhan = 0
for i in list2:
    done = False
    this_cluster_list = []
    length_i = len(i)
    for j in list1:
        intersection = i.intersection(j)
        length_intersection = len(intersection)
        
        if length_intersection == length_i:
            done = True
            break
        
        if length_intersection != 0:
            this_cluster_list.append(j)

    if done == True:
        continue

    tach = True
    print(length_i)
    for j_cluster in this_cluster_list:
        if i.issuperset(j_cluster):
            print(f'subset: {len(j_cluster)}')
        else:
            tach = False
            difference = j_cluster.difference(i)
            intersection = i.intersection(j_cluster)
            print(f"giong: {len(intersection)}")
            print(f"khac: {len(difference)}")
    
    if tach == True:
        chap_nhan += 1
    else:
        khong_chap_nhan +=1

print(chap_nhan)
print(khong_chap_nhan)

print('-------------------------------------------')


ok = 0
no_ok = 0
for j in list1:
    this_cluster_list = []
    for i in list2:
        intersection = i.intersection(j)
        if len(intersection) < len(j) and len(intersection) != 0:
            this_cluster_list.append(i)
    
    if len(this_cluster_list) == 0:
        continue

    tach = True
    print(len(j))
    for i_cluster in this_cluster_list:
        if j.issuperset(i_cluster):
            print(f'subset: {len(i_cluster)}')
        else:
            tach = False
            difference = i_cluster.difference(j)
            intersection = j.intersection(i_cluster)
            print(f"giong: {len(intersection)}")
            print(f"khac: {len(difference)}")

    if tach == True:
        ok += 1
    else:
        no_ok +=1

print(ok)
print(no_ok)