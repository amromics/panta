# chọn random các mẫu vào subset

import os
from glob import glob
import random
import shutil

os.chdir('/home/ntanh1999/amromics/amromics/pan-genome')

gff_files = glob('data/Sp200/*.gff')


# test random set
for j in [50, 100]:
    random.seed(62)
    gff_files.sort()
    for i in range(0, 3):
        
        cmd = f'/usr/bin/time -v python3 pan-genome.py main -c tests/Sp200/test7/{str(j)} --dont-split -t 8 '
        cmd += ' '.join(gff_files[:j])
        os.system(cmd)

        cmd = f'/usr/bin/time -v python3 pan-genome.py add -c tests/Sp200/test7/{str(j)} --dont-split -t 8 '
        cmd += ' '.join(gff_files[j:])
        os.system(cmd)

        cmd = f'python3 experiments/compare_result.py tests/Sp200/test7/{str(j)} tests/Sp200/roary_nosplit'
        os.system(cmd)

        cmd = f'python3 experiments/compare_result_2.py tests/Sp200/test7/{str(j)} tests/Sp200/test7/200'
        os.system(cmd)

        random.shuffle(gff_files)

# # nosplit 
# for j in [50]:
#     random.seed(62)
#     gff_files.sort()
        
#     cmd = f'/usr/bin/time -v python3 pan-genome.py main -c tests/Sp200/test11/{str(j)} -t 8 '
#     cmd += ' '.join(gff_files[:j])
#     os.system(cmd)

#     cmd = f'/usr/bin/time -v python3 pan-genome.py add -c tests/Sp200/test11/{str(j)} -t 8 '
#     cmd += ' '.join(gff_files[j:])
#     os.system(cmd)

#     cmd = f'python3 experiments/compare_result.py tests/Sp200/test11/{str(j)} tests/Sp200/roary_split'
#     os.system(cmd)

#     cmd = f'python3 experiments/compare_result_2.py tests/Sp200/test8/{str(j)} tests/Sp200/200'
#     os.system(cmd)