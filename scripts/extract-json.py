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
def export_json(inp_dir,collectionID,exp_dir):
    updateCollectionHistory(exp_dir,args.id,args.id,"Not Ready") 
    #look for dump file:
    if (os.path.isfile(inp_dir+"/"+collectionID+"_dump.json")):
        print("Dump file not found!")
        return
    dumpfile=os.path.join(inp_dir,+collectionID+"_dump.json")
    exp_dir_current=os.path.join(exp_dir,collectionID)
    if not os.path.exists(exp_dir_current):
        os.makedirs(exp_dir_current)
    report=json.load( open( dumpfile ) )
     
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
    if not os.path.exists(exp_dir_current+"/set"):
        os.makedirs(exp_dir+"/set")
    set_result.append({'group':'phylo_heatmap','data':exportAMRHeatmap(report,exp_dir_current)})
    set_result.append({'group':'pan_sum','data':exportPangenomeSumary(report['set']['pangenome']+'/summary_statistics.txt',exp_dir_current)})
    set_result.append({'group':'pan_cluster','data':exportPangenomeCluster(report['set']['pangenome']+'/gene_presence_absence.csv',exp_dir_current)})
    set_result.append({'group':'phylogeny_tree','data':exportPhylogenyTree(report['set']['phylogeny']+'/parsnp.tree')})
    set_result.append({'group':'gene_alignments','data':exportMultiAlignment(report,exp_dir_current)})
    collection_report={"samples":samples,"results":set_result}
    json.dump( collection_report, open( exp_dir+'/set.json', 'w' ))
    updateCollectionHistory(exp_dir,collectionID,collectionID,"Ready")   

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

def updateCollectionHistory(export_dir,collectionID,collectionName,status):
    if not os.path.exists(export_dir+"/collections.json"):
        #make the empty collection
        emp_arr={collections:[]}
        json.dump(emp_arr,open( export_dir+"/collections.json", 'w' ))
    print(export_dir+"/collections.json")
    collections=json.load(open( export_dir+"/collections.json"))
    isExist=False
    for key in collections["collections"]:
        print(key["collectionID"])
        if key["collectionID"]==collectionID:
            key["collectionName"]=collectionName
            key["status"]=status
            isExist=True
    if not isExist:
        collections["collections"].append({"collectionID":collectionID,"collectionName":collectionName,"status":status})
    json.dump(collections,open( export_dir+"/collections.json", 'w' ))
 
def main(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(
        prog='extract tools for amromics-viz',
        description='extract tools for amromics-viz')
    
    parser.add_argument('--id', help='Colletion ID',type=str)
    parser.add_argument('--inp', help='The output folder of pipeline',type=str)
    parser.add_argument('--out', help='The folder hold exported json data', default='export')

    args = parser.parse_args()
    export_json(args.inp,args.id,args.out)
if __name__ == "__main__":
    main()
