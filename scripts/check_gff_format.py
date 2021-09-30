import argparse
import os
import re
import csv
import logging

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(levelname)s : %(message)s',
    datefmt='%I:%M:%S')
logger = logging.getLogger(__name__)


def rewrite_gff(infile, outfile, ID, fasta=None):
    found_fasta = False
    count = 1
    with open(infile, 'r') as in_fh, open(outfile, 'w') as out_fh:
        for line in in_fh:
            if found_fasta == True:
                out_fh.write(line)
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                out_fh.write(line)
                continue
            if re.match(r"^#", line) != None:
                out_fh.write(line)
                continue
            cells = line.rstrip().split('\t')
            if len(cells) != 9:
                raise Exception(f'{infile} do not follow gff3 format')
            if cells[2] != 'CDS':
                out_fh.write(line)
                continue

            tags = cells[8].split(';')
            found = None
            for i,tag in enumerate(tags):
                result = re.match(r"^ID=(.+)", tag)
                if result != None:
                    found = i
                    break
            if found != None:
                tags[found] = "ID={}_{:05d}".format(ID, count)
                cells[8] = ';'.join(tags)
                out_fh.write('\t'.join(cells) + '\n')
            else:
                tags.append("ID={}_{:05d}".format(ID, count))
                cells[8] = ';'.join(tags)
                out_fh.write('\t'.join(cells) + '\n')
            count += 1
        if found_fasta == False:
            out_fh.write('##FASTA' + '\n')

    
    if found_fasta == False:
        if fasta != None:
            if not os.path.isfile(fasta):
                raise Exception(f'{fasta} does not exist')
            os.system (f'cat {fasta} >> {outfile}')
        else:
            os.remove(outfile)
            raise Exception(f'{infile} requires fasta')

    logging.info(f'Write {outfile}')


def main():
    parser = argparse.ArgumentParser(description='Check and fix gff files')
    parser.add_argument('-o', '--outdir', help='output directory', required=True, type=str)
    parser.add_argument('-t', '--tsv', help='tsv input file',default=None, type=str)
    parser.add_argument('-g', '--gff', help='gff input files',default=None, nargs='*', type=str)
    args = parser.parse_args()

    outdir = args.outdir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if args.tsv != None:
        with open(args.tsv,'r') as fh:
            csv_reader = csv.reader(fh, delimiter='\t')
            for row in csv_reader:
                ID = row[0]
                outfile = os.path.join(outdir, ID + '.gff')
                fasta = row[2]
                if row[2] == '':
                    fasta = None
                rewrite_gff(row[1], outfile, ID, fasta)
    elif args.gff != None:
        gff_list = args.gff
        for gff in gff_list:
            if not gff.endswith('gff'):
                raise Exception(f'{gff} should be a gff3 file')
            base_name = os.path.basename(gff)
            ID = base_name.rsplit('.', 1)[0]
            outfile = os.path.join(outdir, ID + '.gff')
            rewrite_gff(gff, outfile, ID)
    else:
        raise Exception(f'Please specify -t or -g')
    


if __name__ == "__main__":
    main()


    # python3 check_gff.py -o ~/out -t /home/ntanh1999/gb.tsv

    # python3 check_gff.py -o ~/out -g /home/ntanh1999/disk/datasets/data/Sp200/*.gff

    # tsv file contains 3 collums: sample_id; gff_file; fasta_file