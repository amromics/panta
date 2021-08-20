import os

cmd = '/usr/bin/time -v panaroo -i /home/ntanh1999/pan-genome/amromics/amromics/pan-genome/data/Sp25/*.gff -o /home/ntanh1999/pan-genome/panaroo/Sp25 -t 8 --clean-mode strict > /home/ntanh1999/pan-genome/panaroo/Sp25/Sp25.log  2>&1'
os.system(cmd)
cmd = '/usr/bin/time -v panaroo -i /home/ntanh1999/amromics/amromics/pan-genome/data/Sp400/*.gff -o /home/ntanh1999/panaroo/Sp400 -t 8 --clean-mode strict > /home/ntanh1999/panaroo/Sp400/Sp400.log  2>&1'
os.system(cmd)
cmd = '/usr/bin/time -v panaroo -i /home/ntanh1999/pan-genome/amromics/amromics/pan-genome/data/Sp616/*.gff -o /home/ntanh1999/pan-genome/panaroo/Sp616 -t 8 --clean-mode strict > /home/ntanh1999/pan-genome/panaroo/Sp616/Sp616.log  2>&1'
os.system(cmd)