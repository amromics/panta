import json
import os

os.chdir('/home/ntanh1999/amromics/amromics/pan-genome')

gene_position = json.load(open(os.path.join('tests/Kp26/test3/gene_position.json'), 'r'))
inflated_clusters = json.load(open(os.path.join('tests/Kp26/test3/inflated_clusters.json'), 'r'))

gene_to_cluster_index = {gene:index for index, genes in enumerate(inflated_clusters) for gene in genes}

new_gene_position = {}

with open('aaaa', 'w') as fh:
    for sample in gene_position:
        fh.write('#'+sample + '\n')
        new_gene_position[sample] = {}
        for contig in gene_position[sample]:
            fh.write('>'+contig + '\n')
            new_gene_position[sample][contig] = []
            for gene in gene_position[sample][contig]:
                fh.write(str(gene_to_cluster_index[gene])+'_')
                new_gene_position[sample][contig].append(gene_to_cluster_index[gene])
            fh.write('\n')

group = [19, 171,733, 7143]
group_ls = set()
with open('bbbb', 'w') as fh:
    for sample in new_gene_position:
        for contig in new_gene_position[sample]:
            genes_of_contig = new_gene_position[sample][contig]
            indexes = [i for i,x in enumerate(genes_of_contig) if x in group]
            index_ls = []
            for i in indexes:
                pre_index = i - 5
                post_index = i + 6
                if pre_index < 0:
                    pre_index = 0
                length_of_contig = len(genes_of_contig)
                if post_index >= length_of_contig:
                    post_index = length_of_contig
                index_ls.extend(range(pre_index,post_index))
            index_ls.sort()
            ranges =[]
            for x in index_ls:
                if not ranges:
                    ranges.append([x])
                elif x-prev_x == 0:
                    pass
                elif x-prev_x == 1:
                    ranges[-1].append(x)
                else:
                    ranges.append([x])
                prev_x = x
            for x in ranges:
                neighbour_genes = genes_of_contig[x[0]:x[-1]+1]
                for g in neighbour_genes:
                    group_ls.add(g)
                    fh.write(str(g)+'_')
                fh.write('\t'+contig)
                fh.write('\n')
    