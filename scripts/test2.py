# thêm mẫu nhiều lần

import os
from glob import glob
import random
import shutil

os.chdir('/home/ntanh1999/amromics/amromics/pan-genome')

gff_files = glob('data/Sp200/*.gff')


# test for add sample many times
random.seed(62)
gff_files.sort()
for i in range(0, 10):
    random.shuffle(gff_files)

    cmd = '/usr/bin/time -v python3 pan-genome.py main -c tests/Sp200/test4/10010 --dont-split -t 8 '
    cmd += ' '.join(gff_files[:100])
    os.system(cmd)

    for i in range(100,200,10):
        cmd = '/usr/bin/time -v python3 pan-genome.py add -c tests/Sp200/test4/10010 --dont-split -t 8 '
        cmd += ' '.join(gff_files[i:i+10])
        os.system(cmd)

    cmd = f'python3 experiments/compare_result.py tests/Sp200/test4/10010 tests/Sp200/Sp200_roary_nosplit'
    os.system(cmd)

