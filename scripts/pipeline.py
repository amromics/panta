import sys
import argparse
import logging
import json
import csv
import os, shutil
from Bio import SeqIO
import datetime
import pandas as pd

NUM_CORES_DEFAULT = 8

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

def get_assembly(sample,base_dir):
    """

    :param sample:
    :param base_dir:
    :return:
    """
    path_out = os.path.join(base_dir, sample['id'] + '_assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    contigs = list(SeqIO.parse(sample['files'], "fasta"))
    assembly_file = os.path.join(path_out, sample['id'] + '_contigs.fasta')
    contigs = sorted(contigs, key=len, reverse=True)

    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id = sample['id']+'_C'+str(i)
            contig.description = ''
            SeqIO.write(contig, f, "fasta")
    return assembly_file


def run_single_sample(sample, base_dir='.', threads=0, memory=50, trim=False,timing_log=None):
    sample['execution']['start'] = str(datetime.datetime.now())

    if sample['type'] != 'asm':
        raise Exception('Only support asm type!')
        # sample['execution']['out']['assembly'] = assemble_skesa(sample, base_dir=base_dir, threads=0, memory=memory,timing_log=timing_log)
    else:
        sample['execution']['out']['assembly'] = get_assembly(sample, base_dir=base_dir)

    sample['execution']['out']['annotation'] = annotate_prokka(sample, base_dir=base_dir, timing_log=timing_log, threads=threads)
    sample['execution']['out']['mlst'] = mlst(sample, base_dir=base_dir, threads=threads)
    sample['execution']['out']['resistome'] = detect_amr(sample, base_dir=base_dir, timing_log=timing_log, threads=threads)
    sample['execution']['out']['virulome'] = detect_virulome(sample, base_dir=base_dir, threads=threads)
    sample['execution']['out']['plasmid'] = detect_plasmid(sample, base_dir=base_dir, threads=threads)
    sample['execution']['end'] = str(datetime.datetime.now())
    return sample


def annotate_prokka(sample, base_dir='.', timing_log=None, threads=0, overwrite=False):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'prokka_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    file_out = os.path.join(path_out, sample['id'] + '.gff')
    if os.path.isfile(file_out) and (not overwrite):
        # Dont run again if gff file exists
        print('GFF file found, skip annotating')
        return path_out

    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=sample['id'], path=path_out)
    if sample['genus']:
        cmd += ' --usegenus --genus ' + sample['genus']
    if sample['species']:
        cmd += ' --species ' + sample['species']
    if sample['strain']:
        cmd += ' --strain ' + sample['strain']

    # Disable this for now so that we dont have to install signalp
    # if sample['gram']:
    #    cmd += ' --gram ' + sample['gram']

    cmd += ' ' + sample['execution']['out']['assembly']
    if run_command(cmd, timing_log) != 0:
        raise Exception('Command {} returns non-zero!'.format(cmd))
        return None
    return path_out


def mlst(sample,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'mlst_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    mlst_out = os.path.join(path_out, sample['id'] + '_mlst.tsv')
    cmd = 'mlst --quiet --threads {threads} --nopath {infile} > {outfile}'.format(threads=threads,infile=sample['execution']['out']['assembly'],outfile=mlst_out)
    if run_command(cmd, timing_log) != 0:
        return None
    return mlst_out


def detect_amr(sample,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'abricate_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    # TODO: replace by consensus db later
    amr_out = os.path.join(path_out, sample['id'] + '_resistome.tsv')
    cmd = 'abricate --quiet --threads {threads} --nopath --db card {infile} > {outfile}'.format(threads=threads,infile=sample['execution']['out']['assembly'],outfile=amr_out)
    if run_command(cmd, timing_log) != 0:
        return None
    return amr_out


###Virulome profiling using abricate with VFDB
def detect_virulome(sample,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'abricate_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    vir_out = os.path.join(path_out, sample['id'] + '_virulome.tsv')
    cmd = 'abricate --quiet --threads {threads} --nopath --db vfdb {infile} > {outfile}'.format(threads=threads,infile=sample['execution']['out']['assembly'],outfile=vir_out)
    if run_command(cmd, timing_log) != 0:
        return None
    return vir_out


def detect_plasmid(sample,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'abricate_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    #Plasmid finder
    oriREP_out = os.path.join(path_out, sample['id'] + '_plasmid.tsv')
    cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(threads=threads,infile=sample['execution']['out']['assembly'],outfile=oriREP_out)
    if run_command(cmd, timing_log) != 0:
        return None
    return oriREP_out


def run_roary(report,threads=0, base_dir='.', timing_log=None, overwrite=False):
    """
        Run roay make pangeome analysis (using prokka results in previous step)
        :param report: result holder
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """
    gff_list = []
    for sample_id in report['samples']:
        sample = report['samples'][sample_id]
        gff_file = os.path.join(sample['execution']['out']['annotation'], sample_id + '.gff')
        assert os.path.isfile(gff_file)
        gff_list.append(gff_file)

    roary_folder = os.path.join(base_dir,'set/roary')
    roary_output = os.path.join(roary_folder, 'core_alignment_header.embl')
    if os.path.isfile(roary_output) and (not overwrite):
        print('roary has run')
        return report

    #Make sure the directory is not there or roary will add timestamp
    if os.path.isfile(roary_folder):
        os.remove(roary_folder)
    if os.path.exists(roary_folder):
        shutil.rmtree(roary_folder)

    myCmd = 'roary -p {} -f {} -e -n -v '.format(threads, roary_folder) + ' '.join(gff_list)
    run_command(myCmd, timing_log)
    report['set']['pangenome'] = roary_folder
    return report


def run_phylogeny(report,base_dir,threads=0, timing_log=None):
    """
        Run parsnp to create phylogeny tree
        :param report: result holder
        :param ref_genome: path to reference genome, if equal None, one of genome in genome directory will be choosed to be reference.
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """
    phylogeny_folder=os.path.join(base_dir,'set/phylogeny')
    if not os.path.exists(phylogeny_folder):
        os.makedirs(phylogeny_folder)
    genome_dir=os.path.join(base_dir,'temp/fasta')
    if not os.path.exists(genome_dir):
        os.makedirs(genome_dir)

    for sample in report['samples']:
        shutil.copy(report['samples'][sample]['execution']['out']['assembly'], genome_dir+'/'+os.path.basename(report['samples'][sample]['execution']['out']['assembly'] ))

    # take first genome to get reference genome
    files=os.listdir(genome_dir)
    candidate_ref=None
    for f in files :
        if f.endswith(('.fna', '.fa', '.fn', '.fasta')) :
            #check if seq.desc contain '-' character, may cause error with parsnp
            containSpecialCharacter=False
            for seq in SeqIO.parse(os.path.join(genome_dir,f), "fasta"):
                if '-' in seq.id:
                    containSpecialCharacter=True
                    break
            if containSpecialCharacter:
                # keep checking remain genome to find ref
                continue
            else:
                candidate_ref=os.path.join(genome_dir,f)
                break
    if candidate_ref ==None:
        print('Cannot determine appropriate reference genome, may be description of contigs contain special characters (-)')
    else:
        ref_genome=candidate_ref
    myCmd = 'parsnp -p {} -d {} -r {} -o {}'.format(threads, genome_dir, ref_genome, phylogeny_folder)
    run_command(myCmd, timing_log)

    report['set']['phylogeny'] = phylogeny_folder
    return report


def runAlignment(report,base_dir, timing_log=None):
    gene_cluster_file=report['set']['pangenome']+'/gene_presence_absence.csv'
    dict_cds={}
    for id in report['samples']:
        for seq in SeqIO.parse(report['samples'][id]['execution']['out']['annotation']+'/'+id+'.ffn', "fasta"):
            dict_cds[seq.id]=seq

    #make folder contains sequences for each gene
    alignment_dir=os.path.join(base_dir,'set/alignments')
    n_col=0
    fieldnames=[]
    f = open(gene_cluster_file,'r')
    reader = csv.reader(f, delimiter=',')
    fieldnames=next(reader)
    
    for row in reader:
        #extract
        gene_dir=os.path.join(base_dir, 'temp/alignments/genes/'+row[0]) 
        gene_file=''
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)
        for i in range(len(fieldnames)):
            if i>13:               
                if(not row[i]==''):
                    gene_file = os.path.join(gene_dir, fieldnames[i] + '.fasta')
                    
                    SeqIO.write(dict_cds[row[i].split('\t')[0]],gene_file,"fasta")
        if not os.path.exists(os.path.join(alignment_dir, row[0])):
            os.makedirs(os.path.join(alignment_dir, row[0]))
        myCmd = 'parsnp -p {} -d {} -r {} -o {}'.format(2,gene_dir,gene_file,os.path.join(base_dir, 'set/alignments/'+row[0]))
        run_command(myCmd, timing_log)
    f.close()
    report['set']['alignments']=alignment_dir
    return report


def pipeline_func(args):
    report = {'samples': {}, 'set': {}}
    workdir=args.work_dir+"/"+args.id
    sample_df = pd.read_csv(args.input, sep='\t')
    sample_df.fillna('', inplace=True)
    for _, row in sample_df.iterrows():
        sample = {}
        sample['id'] = str(row['Sample ID'])
        sample['name'] = row['Sample Name']
        sample['type'] = row['Input Type']
        sample['files'] = row['Files']
        sample['genus'] = row['Genus']
        sample['species'] = row['Species']
        sample['strain'] = row['Strain']
        sample['gram'] = row['Gram']
        metadata = row['Metadata'].split(';')
        mt = {}
        if len(metadata) > 0:
            for kv in metadata:
                if len(kv.split(':')) == 2:
                    k, v = kv.split(':')
                    mt[k] = v
        sample['metadata'] = mt
        sample['execution'] = {}
        report['samples'][sample['id']] = sample

    #run single sample pipeline
    for id in report['samples']:
        sample_dir=workdir+'/samples/'+str(id)
        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)
        report['samples'][id]['execution']['out']={}
        report['samples'][id] = run_single_sample(report['samples'][id],base_dir=sample_dir, threads=args.threads, memory=args.memory, timing_log=args.time_log)
    
    report=run_roary(report,threads=args.threads,base_dir=workdir, timing_log=args.time_log)
    report=run_phylogeny(report,threads=args.threads,base_dir=workdir, timing_log=args.time_log)
    report=runAlignment(report,base_dir=workdir, timing_log=args.time_log)
    json.dump(report, open( workdir+"/"+args.id+"_dump.json", 'w' ))
    #clean up
    if os.path.exists(workdir+"/temp"):
        shutil.rmtree(workdir+"/temp")


def version_func():
    print('V.1.0')


def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='pipeline_bacterial_analysis',
        description='Tool for amrb pipeline')
    subparsers = parser.add_subparsers(title='sub command', help='sub command help')

    version_cmd = subparsers.add_parser(
        'version', description='Print version of this and other binaries',
        help='Print version of this and other binaries',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    version_cmd.set_defaults(func=version_func)

    pa_cmd = subparsers.add_parser(
        'pa', description='NGS analysis pipeline', help='NGS analysis pipeline',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa_cmd.set_defaults(func=pipeline_func)
    pa_cmd.add_argument('--id', help='Colletion ID',type=str)
    pa_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    pa_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    pa_cmd.add_argument('-i', '--input', help='Input file',type=str)
    pa_cmd.add_argument('--work-dir', help='Working folder', default='out')
    pa_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)

    args = parser.parse_args(arguments)
    return args.func(args)

if __name__ == "__main__":
    main()
