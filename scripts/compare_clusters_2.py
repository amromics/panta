# compare 2 pipeline
# compare base on clusters_2


import csv
import sys
file1 = sys.argv[1]
file2 = sys.argv[2]

size_collection = int(sys.argv[3])

def get_clusters(result_file):
    clusters = []
    count_samples = []
    reader = csv.reader(open(result_file, 'r'))
    next(reader)
    for row in reader:
        this_cluster = set()
        count = 0
        for cell in row[8:]:
            if cell != "":
                count += 1

            genes = cell.split('\t')
            for gene in genes:
                if gene == '':
                    continue
                this_cluster.add(gene)
        clusters.append(this_cluster)
        count_samples.append(count)
    return clusters, count_samples


clusters_1, count_1 = get_clusters(file1)
clusters_2, count_2 = get_clusters(file2)

dictionary ={'core':{'same':0, 'diff':0}, 'soft':{'same':0, 'diff':0}, 'shell': {'same':0, 'diff':0}, 'cloud': {'same':0, 'diff':0}}
total = 0
for i, count in zip(clusters_2,count_2):
    found = False
    percent = count / size_collection
    if percent >= 0.99:
        key = 'core'
    elif percent >= 0.95:
        key = 'soft'
    elif percent >= 0.15:
        key = 'shell'
    else:
        key = 'cloud'

    for j in clusters_1:
        if i == j:
            dictionary[key]['same'] += 1
            found = True
            break
    
    if found == False:
        dictionary[key]['diff'] += 1
    
    total += 1

print(dictionary)
print(total)