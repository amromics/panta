import sys
import argparse
import logging
import json
import subprocess
import csv
import multiprocessing
import os, shutil, glob
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import datetime
import base64
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
    #cmd = "bash -c '{}'".format(cmd)
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
    #cmd = "bash -c '{}'".format(cmd)
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
    #cmd = "bash -c '{}'".format(cmd)
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
    #cmd = "bash -c '{}'".format(cmd)
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
    #cmd = "bash -c '{}'".format(cmd)
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

    # gff_folder=os.path.join(base_dir,'temp/gff')
    # if not os.path.exists(gff_folder):
    #     os.makedirs(gff_folder)
    # # assumtion that the path is : base_dir/sample_id/output
    # for sample in report['samples']:
    #     for root, dirs, files in os.walk(report['samples'][sample]['execution']['out']['annotation']):
    #         for _file in files:
    #             if _file.endswith('.gff'):
    #                 #print ('Found file in: ' , str(root),',',_file)
    #                 newname = sample+'.gff'
    #                 #print ('new  file: ' , newname)
    #                 shutil.copy(os.path.abspath(str(root) + '/' + _file), gff_folder+'/'+newname)

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
    #take first genome to get reference genome
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
                    #keep checking remain genome to find ref
                continue
            else:
                candidate_ref=os.path.join(genome_dir,f)
                break
    if candidate_ref ==None:
        print('Can not determine appropriate reference genome, may be description of contigs contain special characters (-)')
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
    #report=json.load( open( "temp.json" ) )
    
    #json.dump(report, open( "temp.json", 'w' ))
    
    report=run_roary(report,threads=args.threads,base_dir=workdir, timing_log=args.time_log)
    report=run_phylogeny(report,threads=args.threads,base_dir=workdir, timing_log=args.time_log)
    report=runAlignment(report,base_dir=workdir, timing_log=args.time_log)
    json.dump(report, open( "temp.json", 'w' ))
    #report=json.load( open( "temp.json" ))
    export_json(report,args.export_dir)

def export_json(report,exp_dir):

    #export single samples
    samples=[]
    for id in report['samples']:
        out=report['samples'][id]['execution']['out']
        result=[]
        #handle assembly results
        ret_asm={'group':'CONTIG'}
        ret_asm['data']=exportAssembly(out['assembly'])
        result.append(ret_asm)
        #handle mgfglst results   gf
        ret_mlst={'group':'MLST'}        
        ret_mlst['data']=extract_mlst(out['mlst'])
        result.append(ret_mlst)
        ret_vir={'group':'VIR'}
        ret_vir['data']=find_virulome(out['virulome'])
        result.append(ret_vir)
        ret_amr={'group':'AMR'}
        ret_amr['data']=find_amr(out['resistome'])
        result.append(ret_amr)
        ret_plasmid={'group':'PLASMID'}
        ret_plasmid['data']=find_plasmid(out['plasmid'])
        result.append(ret_plasmid)
        ret_annotation={'group':'ANNOTATION'}
        ret_annotation['data']=exportKnowgene(out['annotation'])
        result.append(ret_annotation)
        report['samples'][id]['execution']['result']=result
        save_sample_result(report['samples'][id],exp_dir)
        samples.append({
            'id':report['samples'][id]['id'],\
            'name':report['samples'][id]['name'],
            'type':report['samples'][id]['type'],
            'files':report['samples'][id]['files'], 
            'genus':report['samples'][id]['genus'],
            'species':report['samples'][id]['species'],
            'strain':report['samples'][id]['strain'],
            'gram':report['samples'][id]['gram'],
            'metadata':report['samples'][id]['metadata']

        })

    set_result=[]
    if not os.path.exists(exp_dir+"/set"):
        os.makedirs(exp_dir+"/set")
    set_result.append({'group':'phylo_heatmap','data':exportAMRHeatmap(report,exp_dir)})
    set_result.append({'group':'pan_sum','data':exportPangenomeSumary(report['set']['pangenome']+'/summary_statistics.txt',exp_dir)})
    set_result.append({'group':'pan_cluster','data':exportPangenomeCluster(report['set']['pangenome']+'/gene_presence_absence.csv',exp_dir)})
    set_result.append({'group':'phylogeny_tree','data':exportPhylogenyTree(report['set']['phylogeny']+'/parsnp.tree')})
    set_result.append({'group':'gene_alignments','data':exportMultiAlignment(report,exp_dir)})
    collection_report={"samples":samples,"results":set_result}
    json.dump( collection_report, open( exp_dir+'/set.json', 'w' ))

def exportAssembly(contigs_file_contents):
     #contigs_file: contigs.fasta
    contigs_stat={}
    contigs_stat['contigs']=[]
    genome_length=0
    pattern='GCgc'
    num_gc=0
    
    seq_dict={}
    skew_list=[]
    content_list=[]
    for seq in SeqIO.parse(contigs_file_contents, "fasta"):
        current_contig=''
        nodename=seq.id
        length=len(seq.seq)
        contigs_stat['contigs'].append({'name':nodename,'length':length})
        current_contig=nodename
        seq_dict[current_contig]=str(seq.seq)
        genome_length=genome_length+length
        for i in range(len(seq.seq)):
            if seq.seq[i] in pattern:
                num_gc=num_gc+1   
        
    #running window 1000 step 10 for skew
    for c in seq_dict.keys():
        
        list_GCskew=[]
        list_GCcontent=[]
        window=1000
        slide=100
        if len(seq_dict[c])<1000:
            window=len(seq_dict[c])
        for i in range(0,len(seq_dict[c])-window+1,slide):
            Seq=seq_dict[c][i:i+window].upper()
            G=Seq.count('G')
            C=Seq.count('C')
            A=Seq.count('A')
            T=Seq.count('T')
            GC=0
            if  G!=0 or C!=0:
                GC=(C-G)/(G+C+0.0)
          
           
            list_GCskew.append(round(GC,3))
            list_GCcontent.append(round((G+C)/(G+C+A+T),3))
        skew_list.append({'contig':c,'GC':list_GCskew})
        content_list.append({'contig':c,'GC':list_GCcontent})


    #call n50          
    s=0
    n50=0
    for i in range(len(contigs_stat['contigs'])):
        s=s+ contigs_stat['contigs'][i]['length']
        if s>=genome_length/2:
            n50=contigs_stat['contigs'][i]['length']
            break
    contigs_stat['n_contig']=len(contigs_stat['contigs'])        
    contigs_stat['genome_length']=genome_length
    contigs_stat['N50']=n50
    contigs_stat['min_length']=contigs_stat['contigs'][len(contigs_stat['contigs'])-1]['length']
    contigs_stat['max_length']=contigs_stat['contigs'][0]['length']
    contigs_stat['GC']=(num_gc/genome_length)*100
    contigs_stat['skew']=skew_list
    contigs_stat['content']={'window':1000,'step':100,'array':content_list}
    return contigs_stat
def extract_mlst(mlst_file):
    """
    Extract mlst data::: from Quang
    """
    ret={}
    with open(mlst_file) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            ret['st']=row[1]+' ST['+row[2]+']'
            ret['hits']=[]

            for i in range(3,len(row)):
                ret['hits'].append({'locus':row[i].split('(')[0],'allele':row[i].split('(')[1].replace(')','')})
    return ret
def find_virulome(virulome_file):
    set_vir=set()
    ret={}
    ret['hits']=[]
    with open(virulome_file) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            vir={}
            if '~~~' in row['GENE']:
                gene=row['GENE'].split('~~~')[1]

            else:
                gene=row['GENE'].strip()
            set_vir.add(gene)
            vir['sequence']=row['SEQUENCE']
            vir['start']=row['START']
            vir['end']=row['END']
            vir['strain']=row['STRAND']
            vir['gene']=gene
            vir['coverage']=row['%COVERAGE']+'%'
            vir['identity']=row['%IDENTITY']+'%'
            vir['db']=row['DATABASE']
            vir['accession']=row['ACCESSION']
            vir['product']=row['PRODUCT']
            ret['hits'].append(vir)

    s_gene=''

    for v in set_vir:
        s_gene=s_gene+v+', '
    ret['genes']=s_gene[:-2]

    return ret
def find_amr(amr_file):
    set_amr=set()
    ret={}
    ret['hits']=[]
    with open(amr_file) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            amr={}
            set_amr.add(row['RESISTANCE'].strip())
            amr['sequence']=row['SEQUENCE']
            amr['start']=row['START']
            amr['end']=row['END']
            amr['strain']=row['STRAND']
            amr['gene']=row['GENE'].strip()
            amr['coverage']=row['%COVERAGE']+'%'
            amr['identity']=row['%IDENTITY']+'%'
            amr['db']=row['DATABASE']
            amr['accession']=row['ACCESSION']
            amr['product']=row['PRODUCT']
            amr['resistance']=row['RESISTANCE']
            ret['hits'].append(amr)
    str_gene=''

    for v in set_amr:
        str_gene=str_gene+v+', '
    ret['antibiotics']=str_gene[:-2]
    return ret

def find_plasmid(plasmid_file):
    set_plasmid=set()
    with open(plasmid_file) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            set_plasmid.add(row['GENE'].strip())
    ret=''

    for v in set_plasmid:
        ret=ret+v+', '
    ret=ret[:-2]
    return ret
def exportKnowgene(annotation_folder):
    #look for gff file
    gff_file=None
    for root, dirs, files in os.walk(annotation_folder):
        for _file in files:
            if _file.endswith('.gff'):
                gff_file=os.path.abspath(str(root) + '/' + _file)
    if gff_file==None:
        return ''
    knowgene={}
    knowgene['genes']=[]
    f = open(gff_file)
    line = f.readline()
    while line:
        if not line.startswith('##'):
            token=line.split('\t')
            if token[1]=='prokka':
                #start new gene record
                #collect position,strain and gene name:
                contig=token[0]
                start=int(token[3])
                end=int(token[4])
                strain=token[6]
                tok_des=token[8].split(';')
                name=''
                for s in tok_des:
                    if s.startswith('Name='):
                        name=s.split('=')[1]
                #collect type and product
                line = f.readline()
                token2=line.split('\t')
                gene_type=token2[2]
                tok_des2=token2[8].split(';')
                product=''
                for s in tok_des2:
                    if s.startswith('product='):
                        product=s.split('=')[1].strip().replace('\'','')
                hit={}
                hit['contig']=contig
                hit['start']=start
                hit['end']=end
                hit['strain']=strain
                hit['name']=name
                hit['type']=gene_type
                hit['product']=product
                knowgene['genes'].append(hit)
        if line.startswith('##FASTA'):
              break
        #next line
        line = f.readline()
    f.close()
    return knowgene
def save_sample_result(sample,exp_dir):
    sample_dir=os.path.join(exp_dir,'samples')
    if not os.path.exists(sample_dir):
        os.makedirs(sample_dir)
    json.dump( sample, open( sample_dir+'/'+sample['id']+".json", 'w' ))

def exportPangenomeSumary(summary_file,exp_dir):
    f = open(summary_file)
    ret={}
    ret['group']=[]
    lines=f.readlines()
    for i in range(len(lines)-1) :
        tok=lines[i].strip().split("\t")
        ret['group'].append({'name':tok[0],'des':tok[1],'num':tok[2]})    
    ret['total']=lines[len(lines)-1].strip().split('\t')[2]
    f.close()
    save_path=exp_dir+"/set/pangenome_summary.json"
    json.dump( ret, open( save_path, 'w' ))
    
    return "/set/pangenome_summary.json"
def exportPangenomeCluster(pre_abs_file,exp_dir):
    ret={}
    ret['genes']=[]
    with open(pre_abs_file) as tsvfile:
        reader = csv.DictReader(tsvfile,delimiter=',', dialect='excel-tab')
        for row in reader:
            gene={}           
            gene['gene']=row['Gene']
            gene['annotation']=row['Annotation']
            gene['noisolates']=row['No. isolates']
            gene['nosequences']=row['No. sequences']         
            gene['length']=row['Avg group size nuc']
            ret['genes'].append(gene)
    save_path=exp_dir+"/set/pangenome_cluster.json"
    json.dump( ret, open( save_path, 'w' ))
    
    return "/set/pangenome_cluster.json"
    #return ret
def exportAMRHeatmap(report,exp_dir):
    set_genes=set()
    set_class=set()
    set_samples=set()
    heatmap_stats={}
    heatmap_stats['hits']=[]
    set_vir_gene=set()
    for id in report['samples']:
        ret=find_amr(report['samples'][id]['execution']['out']['resistome'])
        for v in ret['hits']:
            set_genes.add(v['gene'])
            set_class.add(v['resistance'])
            hit={}
            hit['sample']=id
            hit['gene']=v['gene']
            hit['type']='amr'
            hit['class']=v['resistance']
            hit['identity']=float(v['identity'].replace('%',''))
            hit['product']=v['product']
            heatmap_stats['hits'].append(hit)
        
        ret=find_virulome(report['samples'][id]['execution']['out']['virulome'])
        for v in ret['hits']:
            set_vir_gene.add(v['gene'])
            #set_class.add(v['resistance'])
            hit={}
            hit['sample']=id
            hit['gene']=v['gene']
            hit['type']='vir'
            #hit['class']=v['class']
            hit['identity']=float(v['identity'].replace('%',''))
            hit['product']=v['product']
            heatmap_stats['hits'].append(hit)
    heatmap_stats['list_genes']=list(set_genes.union(set_vir_gene))
    heatmap_stats['list_class']=list(set_class)
    heatmap_stats['samples']=list(set_samples)
    save_path=exp_dir+"/set/amrheatmap.json"
    json.dump( heatmap_stats, open( save_path, 'w' ))
    return "/set/amrheatmap.json"
def exportPhylogenyTree(treefile):
    data=''
    with open(treefile) as myfile:
        data="".join(line.rstrip() for line in myfile)
    message_bytes = data.encode('ascii')
    base64_bytes = base64.b64encode(message_bytes)
    base64_message = base64_bytes.decode('ascii')
    return base64_message
def exportMultiAlignment(report,exp_dir):
    list_genes=os.listdir( report['set']['alignments'] )
    alignments={}
    alignments['alignments']=[]
    for gene in list_genes:
        if os.path.isdir(report['set']['alignments']+'/'+gene):
            if not os.path.isfile(report['set']['alignments']+'/'+gene+'/parsnp.tree'): 
                continue
            tree=exportPhylogenyTree(report['set']['alignments']+'/'+gene+'/parsnp.tree')
            aln={'gene':gene}
            aln['tree']=tree
            aln['samples']=exportAlignment(gene,report['set']['alignments']+'/'+gene+'/parsnp.xmfa',exp_dir,)
            alignments['alignments'].append(aln)

    return alignments
def exportAlignment(gene,file_xmfa,exp_dir):
    f = open(file_xmfa)
    aligments=[]
    s_dict={}
    recent_index=0
    current_index=0
    seq=''
    line = f.readline()
    while line:
        if line.startswith('#'):
            
            if line.startswith('##SequenceIndex'):
                t=line.strip().split(' ')
                if not t[1] in s_dict:
                    s_dict[t[1]]={}
                    recent_index=t[1]
            if line.startswith('##SequenceFile'): 
                t=line.strip().split(' ')
                sampleid=t[1].replace('.fasta.ref','')
                sampleid=sampleid.replace('.fasta','')
                s_dict[recent_index]['id']=sampleid
        elif line.startswith('>'):
            t=line.split(' ')
            current_index=t[0].split(':')[0].replace('>','')
            
            seq=''
        else:
            seq=seq+line.strip()
            
            s_dict[current_index]['seq']=seq
        #next line
        line = f.readline()
    f.close()
    for sid in s_dict:
        sample={'sample':s_dict[sid]['id'],'seq':s_dict[sid]['seq'].upper().replace('=','')}
        aligments.append(sample)
    if not os.path.exists(exp_dir+"/set/alignments/"):
        os.makedirs(exp_dir+"/set/alignments/")
    save_path=exp_dir+"/set/alignments/"+gene+".json"
    json.dump( aligments, open( save_path, 'w' ))
    
    return "/set/alignments/"+gene+".json"
    #return aligments


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
    pa_cmd.add_argument('-e', '--export-dir', help='Export folder', default='export')
    pa_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)

    args = parser.parse_args(arguments)
    return args.func(args)
if __name__ == "__main__":
    main()
