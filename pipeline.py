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
def run_command(cmd, timing_log=None):
    """
    Run a command line, return the returning code of the command
    :param cmd:
    :param timing_log:
    :return:
    """
    if timing_log is not None:
        cmd = '/usr/bin/time --append -v -o {} bash -c "{}"'.format(timing_log, cmd)
    #logger.info('Running "{}'.format(cmd))
    print(cmd)
    ret = os.system(cmd)
    return ret

###Trimming using trimmomatic
def trim_pe_trimmomatic(report, threads=0, base_dir='.', timing_log=None, **kargs):
    """
    report is a dictionary with field `sample_id`
    :param report:
    :param threads:
    :param base_dir:
    :param kargs:
    :return:
    """
    if threads <= 0:
        threads = NUM_CORES_DEFAULT

    out_dir = os.path.join(base_dir, 'trimmomatic_pe')
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    out_p1 = os.path.join(out_dir, report['sample_id'] + '_R1.fastq.gz')
    out_p2 = os.path.join(out_dir, report['sample_id'] + '_R2.fastq.gz')

    out_s1 = os.path.join(out_dir, report['sample_id'] + '_S1.fastq')
    out_s2 = os.path.join(out_dir, report['sample_id'] + '_S2.fastq')

    cmd = 'trimmomatic.sh PE -threads {threads}'.format(threads=threads)
    cmd += ' {in_p1} {in_p2} {out_p1} {out_s1} {out_p2} {out_s2}'.format(
        in_p1=report['pe1'], in_p2=report['pe2'],
        out_p1=out_p1, out_s1=out_s1, out_p2=out_p2, out_s2=out_s2)
    ret = run_command(cmd, timing_log)
    # Combine single-ended reads into one
    out_s = os.path.join(out_dir, report['sample_id'] + '_S.fastq')
    with open(out_s, 'w') as fn:
        for seq in bioseq.read_sequence_file(out_s1):
            fn.write(seq.format_fastq())
        for seq in bioseq.read_sequence_file(out_s2):
            fn.write(seq.format_fastq())

    if ret == 0:
        report['pe1'] = out_p1
        report['pe2'] = out_p2
        report['se'] = out_s
        return report

###NGS assembly using SPAdes
def assemble_spades(sample, base_dir = '.', threads=0, memory=50, timing_log=None, **kargs):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    sample_id = sample['id']
    path_out = os.path.join(base_dir, sample_id + '_spades')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    cmd = 'spades.py -m {memory} -t {threads} -k 77,99,127 --careful -o {path_out}'.format(
        memory=int(memory), threads=threads, path_out=path_out)
    if len(ret_sample['files'].split(';'))>0 :
        pe1=ret_sample['files'].split(';')[0]
        pe2=ret_sample['files'].split(';')[1]
        cmd += ' --pe1-1 {pe1} --pe1-2 {pe2}'.format(pe1=pe1, pe2=pe2)
    else:
        cmd += ' --s1 {se}'.format(se=sample['files'])

    ret = run_command(cmd, timing_log)
    if ret != 0:
        return None

    #Read in list of contigs
    contigs = list(SeqIO.parse(os.path.join(path_out, 'contigs.fasta'), "fasta"))
    #TODO: filter based on coverage

    contigs = sorted(contigs, key=len, reverse=True)

    assembly_file = os.path.join(path_out, sample_id + '_contigs.fasta')

    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id=ret_sample['id']+'_C'+str(i)

            SeqIO.write(contig,assembly_file,"fasta")

    return assembly_file
    # TODO: Get the graph
    # Export report:


def assemble_skesa(sample, base_dir = '.', threads=0, memory=50, timing_log=None, **kargs):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    sample_id = sample['id']
    path_out = os.path.join(base_dir, sample_id + '_skesa')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    cmd = 'skesa --memory {memory} --cores {threads} --fastq '.format(
        memory=int(memory), threads=threads)
    if len(ret_sample['files'].split(';'))>0 :
        pe1=ret_sample['files'].split(';')[0]
        pe2=ret_sample['files'].split(';')[1]
        cmd += '{pe1} {pe2}'.format(pe1=pe1, pe2=pe2)
    if 'se' in read_data:
        cmd += '{se}'.format(se=sample['files'])
    assembly_file_raw= os.path.join(path_out, 'contigs.fasta')
    cmd+=' >{output}'.format(output=assembly_file_raw)
    ret = run_command(cmd, timing_log)
    if ret != 0:
        return None

    #Read in list of contigs
    contigs = list(SeqIO.parse(assembly_file_raw, "fasta"))
    #TODO: filter based on coverage
    assembly_file = os.path.join(path_out, sample_id + '_contigs.fasta')
    contigs = sorted(contigs, key=len, reverse=True)
    #logger.info("Read in {} contigs".format(len(contigs)))

    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id=ret_sample['id']+'_C'+str(i)

            SeqIO.write(contig,assembly_file,"fasta")
    #ret_sample['assembly'] = assembly_file
    return assembly_file
def get_assembly(sample,base_dir):
    path_out = os.path.join(base_dir, sample['id'] + '_assembly')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    contigs = list(SeqIO.parse(sample['files'], "fasta"))
    #TODO: filter based on coverage
    print(len(contigs))
    assembly_file = os.path.join(path_out, sample['id'] + '_contigs.fasta')
    contigs = sorted(contigs, key=len, reverse=True)
    #logger.info("Read in {} contigs".format(len(contigs)))

    with open(assembly_file, 'w') as f:
        for i, contig in enumerate(contigs):
            contig.id=sample['id']+'_C'+str(i)
            #print(len(contig))
            SeqIO.write(contig,f,"fasta")
    #ret_sample['assembly'] = assembly_file
    return assembly_file
###Annotation using prokka
def annotate_prokka( sample, base_dir='.', timing_log=None, threads=0):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'prokka_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    cmd = 'prokka --force --cpus {threads} --addgenes --mincontiglen 200'.format(threads=threads)
    cmd += ' --prefix {sample_id} --locus {sample_id} --outdir {path} '.format(sample_id=sample['id'], path=path_out)
    if not sample['genus']=='':
        cmd += ' --genus ' + sample['genus']
    if not sample['species']=='':
        cmd += ' --species ' + sample['genus']
    if not sample['strain']=='':
        cmd += ' --strain ' + sample['strain']
    if not sample['gram']=='':
        cmd += ' --gram ' + sample['gram']
    cmd += ' ' + sample['execution']['out']['assembly']
    cmd = "bash -c '{}'".format(cmd)
    if run_command(cmd, timing_log) != 0:
        return None
    #ret_sample['annotation'] = path_out
    return path_out

###Sequence typing using mlst
def mlst(sample,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'mlst_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    mlst_out = os.path.join(path_out, sample['id'] + '_mlst.tsv')
    cmd = 'mlst --quiet --threads {threads} --nopath {infile} > {outfile}'.format(threads=threads,infile=sample['execution']['out']['assembly'],outfile=mlst_out)
    cmd = "bash -c '{}'".format(cmd)
    if run_command(cmd, timing_log) != 0:
        return None
    return mlst_out







###AMRomes profiling using abricate
def detect_amr(sample,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'abricate_' + sample['id'])
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    #AMR profiling with CARD. TODO: replace by consensus db later
    amr_out = os.path.join(path_out, sample['id'] + '_resistome.tsv')
    cmd = 'abricate --quiet --threads {threads} --nopath --db card {infile} > {outfile}'.format(threads=threads,infile=sample['execution']['out']['assembly'],outfile=amr_out)
    cmd = "bash -c '{}'".format(cmd)
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
    cmd = "bash -c '{}'".format(cmd)
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
    cmd = "bash -c '{}'".format(cmd)
    if run_command(cmd, timing_log) != 0:
        return None
    return oriREP_out




def run_roary(report,threads=0, base_dir='.',):
    """
        Run roay make pangeome analysis (using prokka results in previous step)
        :param report: result holder
        :param base_dir: working directory
        :param threads: number of core CPU
        :return:
    """
    gff_folder=os.path.join(base_dir,'temp/gff')
    if not os.path.exists(gff_folder):
        os.makedirs(gff_folder)
    # assumtion that the path is : base_dir/sample_id/output
    for sample  in report['samples']:
        for root, dirs, files in os.walk(report['samples'][sample]['execution']['out']['annotation']):
            for _file in files:
                if _file.endswith('.gff'):
                    #print ('Found file in: ' , str(root),',',_file)
                    newname = sample+'.gff'
                    #print ('new  file: ' , newname)
                    shutil.copy(os.path.abspath(str(root) + '/' + _file), gff_folder+'/'+newname)

    roary_folder=os.path.join(base_dir,'set/roary')
    #if not os.path.exists(roary_folder):
    #    os.makedirs(roary_folder)
    myCmd = 'roary -p {} -f {} -e -n -v {}/*.gff'.format(threads,roary_folder,gff_folder)
    run_command(myCmd)
    report['set']['pangenome']=roary_folder
    return report
def run_phylogeny(report,base_dir,threads=0):
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
    for sample  in report['samples']:
        shutil.copy(report['samples'][sample]['execution']['out']['assembly'], genome_dir+'/'+os.path.basename(report['samples'][sample]['execution']['out']['assembly'] ))
    #take first genome to get reference genome
    files=os.listdir(genome_dir)
    candidate_ref=None
    for f in files :
        if f.endswith(('.fna', '.fa', '.fn', '.fasta')) :
            #check if seq.desc contain '-' character, may cause error with parsnp
            containSpecialCharacter=False

            for seq in SeqIO.parse(os.path.join(genome_dir,f), "fasta"):
                if '-' in seq.id :
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
    myCmd = 'parsnp -p {} -d {} -r {} -o {}'.format(threads,genome_dir,ref_genome,phylogeny_folder)
    #print(myCmd)
    run_command(myCmd)

    report['set']['phylogeny']=phylogeny_folder
    return report
def runAlignment(report,base_dir):
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
                    SeqIO.write(dict_cds[row[i]],gene_file,"fasta")
        if not os.path.exists(os.path.join(alignment_dir, row[0])):
            os.makedirs(os.path.join(alignment_dir, row[0]))
        myCmd = 'parsnp -p {} -d {} -r {} -o {}'.format(2,gene_dir,gene_file,os.path.join(base_dir, 'set/alignments/'+row[0]))
        run_command(myCmd)
    f.close()
    report['set']['alignments']=alignment_dir
    return report
def run_single_sample(sample, base_dir='.', threads=0, memory=50, trim=False,timing_log=None):
    #handle assembly input, ignore spades and bwa:
    sample['execution']['start']=str(datetime.datetime.now())

    if not sample['type']=='asm':
        #report = assemble_spades(report, base_dir=base_dir, threads=0, memory=memory,timing_log=timing_log)
        sample['execution']['out']['assembly'] = assemble_skesa(sample, base_dir=base_dir, threads=0, memory=memory,timing_log=timing_log)
    else:
        sample['execution']['out']['assembly']=get_assembly(sample, base_dir=base_dir)
    sample['execution']['out']['annotation'] = annotate_prokka(sample, base_dir=base_dir,timing_log=timing_log, threads=threads)
    sample['execution']['out']['mlst'] = mlst(sample, base_dir=base_dir, threads=threads)
    sample['execution']['out']['resistome'] = detect_amr(sample, base_dir=base_dir,timing_log=timing_log, threads=threads)
    sample['execution']['out']['virulome'] = detect_virulome(sample, base_dir=base_dir, threads=threads)
    sample['execution']['out']['plasmid'] = detect_plasmid(sample, base_dir=base_dir, threads=threads)
    sample['execution']['end']=str(datetime.datetime.now())
    return sample
def pipeline_func(args):
    # report = {'samples': {}, 'set':{}}
    # with open(args.input) as tsvfile:
    #     reader = csv.DictReader(tsvfile, dialect='excel-tab')
    #     for row in reader:
    #         sample={}
    #         sample['id']=row['Sample ID']
    #         sample['name']=row['Sample Name']
    #         sample['type']=row['Input Type']
    #         sample['files']=row['Files']
    #         sample['genus']=row['Genus']
    #         sample['species']=row['Species']
    #         sample['strain']=row['Strain']
    #         sample['gram']=row['Gram']
    #         metadata=row['Metadata'].split(';')
    #         mt={}
    #         if len(metadata)>0:
    #             for kv in metadata:
    #                 if len(kv.split(':'))==2:
    #                     k,v=kv.split(':')
    #                     mt[k]=v
    #         sample['metadata']=mt
    #         sample['execution']={}
    #         report['samples'][sample['id']]=sample


    # #run single sample pipeline
    # for id in report['samples']:
    #     sample_dir=args.work_dir+'/samples/'+id
    #     if not os.path.exists(sample_dir):
    #         os.makedirs(sample_dir)
    #     report['samples'][id]['execution']['out']={}
    #     report['samples'][id]=run_single_sample(report['samples'][id],base_dir=sample_dir, threads=args.threads, memory=args.memory, timing_log=args.time_log)

    # report=run_roary(report,threads=args.threads,base_dir=args.work_dir)
    # report=run_phylogeny(report,threads=args.threads,base_dir=args.work_dir)
    # json.dump(report, open( "temp.json", 'w' ))
    
    report=json.load( open( "temp.json" ) )
    #report=runAlignment(report,base_dir=args.work_dir)
    #json.dump(report, open( "temp.json", 'w' ))
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
        #handle mlst results
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
    #save pangenome result

    set_result=[]
    
    set_result.append({'group':'phylo_heatmap','data':exportAMRHeatmap(report)})
    set_result.append({'group':'pan_sum','data':exportPangenomeSumary(report['set']['pangenome']+'/summary_statistics.txt')})
    set_result.append({'group':'pan_cluster','data':exportPangenomeCluster(report['set']['pangenome']+'/gene_presence_absence.csv')})
    set_result.append({'group':'phylogeny_tree','data':exportPhylogenyTree(report['set']['phylogeny']+'/parsnp.tree')})
    set_result.append({'group':'gene_alignments','data':exportMultiAlignment(report)})
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

def exportPangenomeSumary(summary_file):
    f = open(summary_file)
    ret={}
    ret['group']=[]
    lines=f.readlines()
    for i in range(len(lines)-1) :
        tok=lines[i].strip().split("\t")
        ret['group'].append({'name':tok[0],'des':tok[1],'num':tok[2]})    
    ret['total']=lines[len(lines)-1].strip().split('\t')[2]
    f.close()
    return ret
def exportPangenomeCluster(pre_abs_file):
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
    return ret
def exportAMRHeatmap(report):
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
    return heatmap_stats
def exportPhylogenyTree(treefile):
    data=''
    with open(treefile) as myfile:
        data="".join(line.rstrip() for line in myfile)
    message_bytes = data.encode('ascii')
    base64_bytes = base64.b64encode(message_bytes)
    base64_message = base64_bytes.decode('ascii')
    return base64_message
def exportMultiAlignment(report):
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
            aln['samples']=exportAlignment(report['set']['alignments']+'/'+gene+'/parsnp.xmfa')
            alignments['alignments'].append(aln)
    return alignments
def exportAlignment(file_xmfa):
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
    
    return aligments
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
    pa_cmd.add_argument('-t', '--threads', help='Number of threads to use, 0 for all', default=0, type=int)
    pa_cmd.add_argument('-m', '--memory', help='Amount of memory in Gb to use', default=30, type=float)
    pa_cmd.add_argument('-i', '--input', help='Input file',type=str)
    pa_cmd.add_argument('--work-dir', help='Working folder', default='out')
    pa_cmd.add_argument('-e','--export-dir', help='Export folder', default='export')
    pa_cmd.add_argument('--time-log', help='Time log file', default=None, type=str)


    args = parser.parse_args(arguments)
    return args.func(args)
if __name__ == "__main__":
    main()
