# download dataset
# usage: python download.py [Klebsiella.tsv] [column] [output dir]
# column: 2 = assembly (fna.gz) 
# 3 = annotation (gff.gz)  4 = protein (faa.gz)


import os
import sys

ftp_file = sys.argv[1]
column = int(sys.argv[2])
out_dir = sys.argv[3]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

for line in open(ftp_file, 'r'):
    line = line.rstrip()
    cells = line.split('\t')
    acc_id = cells[0]
    ftp_link = cells[column]
    rsync_link = 'rsync:' + ftp_link.split(':')[1]
    out_file = os.path.join(
        out_dir,acc_id + '.' + ftp_link.split('.')[-2] + '.gz')
    if os.path.exists(out_file):
        continue
    cmd = f'rsync --copy-links --times {rsync_link} {out_file}'
    os.system(cmd)




