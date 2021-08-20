import os
import glob

a = os.listdir('Sp616')
with open('cmd.txt','w') as fh:
    for fasta in a:
        sample_id = fasta.rsplit('.',1)[0]
        path = os.path.join('output', sample_id)
        infile = os.path.join('Sp616', fasta)
        cmd = "prokka --force --cpus 8 --addgenes --mincontiglen 200 --prefix {} --locus {} --outdir {} {} ".format(sample_id, sample_id,path,infile)
        fh.write(cmd + '\n')

