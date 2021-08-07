# compare in pair

import os

with open('cmd.txt','w') as fh:
    for i in range(20,220,20):
        for j in range(20,220,20):
            if j < i:
                continue
            cmd = f'python3 compare_result.py Sp200_add_new_{str(i)} Sp200_add_new_{str(j)} && echo {str(i)} {str(j)}'
            fh.write(cmd + '\n')

os.system('parallel -a cmd.txt')