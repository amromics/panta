# compare 2 pipeline
# compare base on clusters_2


import csv
import sys
file1 = sys.argv[1]
file2 = sys.argv[2]

def get_clusters(result_file):
    clusters = []
    for row in csv.reader(open(result_file, 'r')):
        this_cluster = set()
        if row[0] == 'Gene':
            continue
        for cell in row[8:]:
            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                this_cluster.add(gene)
        clusters.append(this_cluster)
    return clusters


clusters_1 = get_clusters(file1)
clusters_2 = get_clusters(file2)


same = 0
diff = 0
total = 0
for i in clusters_2:
    found = False
    
    for j in clusters_1:
        if i == j:
            same += 1
            found = True
            break
    
    if found == False:
        diff += 1
    
    total += 1

print(str(same) + '\t' + str(diff) + '\t'+ str(total))