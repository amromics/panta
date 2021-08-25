import os
import glob

a = os.listdir('disk/datasets/Pseudomonas_800/fna')
with open('cmd.txt','w') as fh:
    for fasta in a:
        sample_id = fasta.rsplit('.',1)[0]
        path = os.path.join('prokka', sample_id)
        infile = os.path.join('disk/datasets/Pseudomonas_800/fna', fasta)
        outfile = os.path.join('prokka', sample_id, sample_id + '.gff')
        if os.path.exists(outfile):
            print('ok')
            continue
        cmd = "prokka --force --cpus 8 --addgenes --mincontiglen 200 --prefix {} --locus {} --outdir {} {} ".format(sample_id, sample_id,path,infile)
        os.system(cmd)
        #fh.write(cmd + '\n')

# parallel --jobs 0 -a cmd.txt