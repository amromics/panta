import os

os.chdir('/home/ntanh1999')

# for filename in os.listdir('/home/ntanh1999/pan-genome-analysis/data/Kp616'):
#     cells = filename.split('.')
#     new_filename = cells[0] + '.gbk'
#     os.system(f'mv /home/ntanh1999/pan-genome-analysis/data/Kp616/{filename} /home/ntanh1999/pan-genome-analysis/data/Kp616/{new_filename}')


cmd = '/usr/bin/time -v ./panX.py -fn ./data/Kp616 -sl Streptococcus_pneumoniae -dmdc -sitr -t 8  1> Kp616.log 2> Kp616.log'