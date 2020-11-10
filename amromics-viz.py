import sys
import argparse
import os, shutil, glob
def run_command(cmd, timing_log=None):
    """
    Run a command line, return the returning code of the command
    :param cmd:
    :param timing_log:
    :return:
    """
    if timing_log is not None:
        cmd = '/usr/bin/time --append -v -o {} bash -c "{}"'.format(timing_log, cmd)
    print(cmd)
    ret = os.system(cmd)
    return ret
def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='extract tools for amromics-viz',
        description='extract tools for amromics-viz')
    
    parser.add_argument('--id', help='Colletion ID',type=str)
    parser.add_argument('--input', help='The tsv file listing samples',type=str)
    parser.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    parser.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    args = parser.parse_args()
    run_command('python scripts/pipeline.py pa --id '+args.id+' -i '+args.input+' -t 0 -m 16 --work-dir data/output')
    run_command('python scripts/extract-json.py --id '+args.id+' --inp data/output --out web-app/static/data' )
    run_command('cd web-app && live-server --port=3000  --entry-file=index.html' )
    
    #check if amromics-vis is compiled
    
  
if __name__ == "__main__":
    main()