

import re
import subprocess


def valid_id(sid):
    """
    Check if string can be a valid id, that is, it can only the following characters:
        - alphanumerical
        - Underscores
        - Dot
    Parameters
    ----------
    sid

    Returns
    -------

    """
    return re.match(r'^[.\w]+$', sid)


def software_version(software_list=None):
    # List of known software and their version
    cmd_versions = {
        # Basic
        'java': 'java -version 2>&1 | head -n 1',
        'python': 'python --version 2>&1',

        # workflow language
        'cromwell': 'cromwell.sh --version 2>&1',

        # Fundamental tools
        'samtools': 'samtools --version  2>&1| head -n 2 | tr "\n" "  "',
        'blast': 'blastn -version 2>&1 | head -n 1',

        # Assemblers
        'spades': 'spades.py -v 2>&1',
        'skesa': 'skesa --version 2>&1 | tail -1',
        'shovill': 'shovill --version 2>&1',

        # Annotations
        'prokka': 'prokka -version 2>&1',
        'mlst': 'mlst --version  2>&1',
        'abricate': 'abricate --version  2>&1|tr "\n" " " && abricate --list|awk \'BEGIN{printf("| Database: ");}NR>1{printf("%s ",$1)}\'',
        'snippy': 'snippy --version 2>&1',

        # Pangenome tools
        'roary': 'roary --version 2>&1 | tail -n 1',
        'parsnp': 'parsnp --version 2>&1 | tail -1',

        # misc
        'trimmomatic': 'trimmomatic -version 2>&1',

        # Others:
        'nodejs':'node --version 2>&1 | tail -1'
        }
    if software_list is None:
        software_list = cmd_versions.keys()

    for sw in software_list:
        if sw in cmd_versions:
            cmd = cmd_versions[sw]
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            (output, err) = p.communicate()
            ret = p.wait()
            if ret == 0:
                print('{:12}= {}'.format(sw, output.decode().strip()))
            else:
                print('{:12}= {}'.format(sw, 'NOT FOUND'))

    for sw in software_list:
        if sw not in cmd_versions:
            print('Cannot check software version for {}'.format(sw))

