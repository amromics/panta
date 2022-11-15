# download dataset
# usage: python download.py [Klebsiella.tsv] [column] [output dir]
# column: 2 = assembly (fna.gz) 3 = annotation (gff.gz)  4 = protein (faa.gz)


import os
import sys
import pandas as pd
import re
import csv
from glob import glob
from datetime import datetime
import shutil
def downloadAndMerge(link,out_dir):
    os.makedirs(out_dir,exist_ok=True)
    os.makedirs('temp',exist_ok=True)
    gff_fna_file=os.path.join(out_dir, acc_id+".gff")
    if os.path.isfile(gff_fna_file):
        return None
    fna_file_link=os.path.join(link,ftp_link.split('/')[-1]+"_genomic.fna.gz")
    downloaded_fna_file = os.path.join('temp',acc_id+ '.fna.gz')
    downloaded_fna_file_unzip = os.path.join('temp',acc_id+ '.fna')
    cmd = f'wget {fna_file_link} -O {downloaded_fna_file}'
    os.system(cmd)
    gff_file_link=os.path.join(link,ftp_link.split('/')[-1]+"_genomic.gff.gz")
    downloaded_gff_file = os.path.join('temp',acc_id+ '.gff.gz')
    cmd = f'wget {gff_file_link} -O {downloaded_gff_file}'
    downloaded_gff_file_unzip = os.path.join('temp',acc_id+ '.gff')
    os.system(cmd)
    if not os.path.exists(downloaded_fna_file) or not os.path.exists(downloaded_gff_file) :
        return None
    cmd=f'gunzip -c {downloaded_fna_file} > {downloaded_fna_file_unzip}'
    print(cmd)
    os.system(cmd)

    cmd=f'gunzip -c {downloaded_gff_file} > {downloaded_gff_file_unzip}'
    print(cmd)
    os.system(cmd)
    if not os.path.isfile(downloaded_fna_file_unzip) or not os.path.isfile(downloaded_gff_file_unzip):
        return None
    gff_fna_file_temp = os.path.join('temp',acc_id+ '_temp.gff')
    gff_fna_file=os.path.join(out_dir, acc_id+".gff")
    cmd=f'bash -c "cat {downloaded_gff_file_unzip}; echo \'##FASTA\' ; cat {downloaded_fna_file_unzip}" > {gff_fna_file_temp}'

    print(cmd)
    os.system(cmd)
    if validGFF(gff_fna_file_temp):
        normalizeGFF(acc_id,gff_fna_file_temp,gff_fna_file)



    if os.path.exists(gff_fna_file) and not validGFF(gff_fna_file):
        os.remove(gff_fna_file)
    shutil.rmtree('temp')
    return gff_fna_file
def validGFF(gff_file):
    with open(gff_file,'r') as in_fh:
        found_fasta=False
        for line in in_fh:
            if found_fasta == True:
                if not line.startswith('>'):
                    if not re.match(r"^[CAGTNSYKcagtnsyk]+$",line):
                        print("invalid fasta line:"+line+" in "+gff_file)
                        return False
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                continue
            if re.match(r"^#", line) != None:
                continue
            line = line.rstrip('\n')
            cells = line.split('\t')
            if len(cells)<8:
                print("invalid line:"+line +" in "+gff_file)
                return False
    return True
def normalizeGFF(acc_id,gff_file_in,gff_file_out):
    if not gff_file_in.endswith('gff'):
        raise Exception(f'{gff} should be a gff3 file')
    base_name = os.path.basename(gff_file_in)
    sample_id = acc_id
    found_fasta = False
    with open(gff_file_in,'r') as in_fh,open(gff_file_out, 'w') as gff_re:
        last_cds=""
        id_number=0
        pre_id=None
        pre_rewrite_id=None
        for line in in_fh:
            if found_fasta == True:
                gff_re.write(line.replace('/','_'))
                continue
            if re.match(r"^##FASTA", line) != None:
                found_fasta = True
                gff_re.write(line)
                continue
            if re.match(r"^#", line) != None:
                gff_re.write(line)
                continue
            line = line.rstrip('\n')
            cells = line.split('\t')
            if len(cells)<8:
                print(line)
                continue
            c_locus=cells[0]+"-"+cells[3]+"-"+cells[4]
            if not c_locus == last_cds:
                last_cds=c_locus
                id_number=id_number+1
            genid=None
            parentid=None
            ##print(line)
            #print(gff_file_in)

            tags=cells[8].split(';')
            for tag in tags:
                ID = re.match(r"^ID=(.+)", tag)
                if ID != None:
                    genid = ID.group(1)
                Parent = re.match(r"^Parent=(.+)", tag)
                if Parent != None:
                    parentid = Parent.group(1)
            newId=None
            newId=sample_id+"-"+cells[0]+"-"+str(id_number)+"-"+genid
            newdes="ID="+newId+";"
            if not parentid==None and parentid==pre_id:
                newdes=newdes+"Parent="+pre_rewrite_id+";"
            for i in range(2,len(tags)):
                newdes=newdes+tags[i]+";"
            newline=""
            for i in range(len(cells)-1):
                newline=newline+cells[i]+"\t"
            newline=newline+newdes+"\n"
            gff_re.write(newline)
            if not genid ==None:
                pre_id=genid
                pre_rewrite_id=newId
    return gff_file_out

ftp_file = sys.argv[1]

out_dir = sys.argv[2]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)
count=0
df = pd.read_csv(ftp_file)
print(len(df))
#drop sample with no refseq
df=df.dropna(subset=['RefSeq FTP'])
print(len(df))
#drop complete sample with small Size
df = df.drop(df[(df["Size(Mb)"] < 5) & (df["Level"] == "Complete")].index)
print(len(df))
#print(len(df[(df["Level"] == "Complete")&(datetime.strptime(df["Release Date"],"%Y-%m-%dT%H:%M:%SZ").year<=2020)]))
input("Press Enter to continue...")
s={}
num_complete_before_2020=0
count_Q1=0
count_Q2=0
count_Q3=0
count_Q4=0
out_dir_complete=os.path.join(out_dir,"complete")
if not os.path.exists(out_dir_complete):
    os.mkdir(out_dir_complete)
for index, row in df.iterrows():
    print(index)
    #print(row)
    acc_id = row["Assembly"]
    ftp_link = row["RefSeq FTP"]
    size=row["Size(Mb)"]
    print(acc_id)
    print(ftp_link)

    if not ftp_link =="":
        date=datetime.strptime(row["Release Date"],"%Y-%m-%dT%H:%M:%SZ")
        print(date.year)
        print(date.month)
        if not date.year in s.keys():
            s[date.year]={}
            s[date.year]["count"]=0
        s[date.year]["count"]=s[date.year]["count"]+1
        if not date.month in s[date.year].keys():
            s[date.year][date.month]=0
        s[date.year][date.month]=s[date.year][date.month]+1
        isComplete2020=False
        if date.year<=2020 and row["Level"]=="Complete":
            if downloadAndMerge(ftp_link,out_dir_complete)is not None:
                num_complete_before_2020=num_complete_before_2020+1
        if not date.year==2021:
            continue
        else:
            dir=os.path.join(out_dir,str(date.year))
            if date.month<4:
                dir=os.path.join(out_dir,"Q1")
                count_Q1=count_Q1+1
            elif date.month<7:
                dir=os.path.join(out_dir,"Q2")
                count_Q2=count_Q2+1
            elif date.month<10:
                dir=os.path.join(out_dir,"Q3")
                count_Q3=count_Q3+1
            else:
                dir=os.path.join(out_dir,"Q4")
                count_Q4=count_Q4+1
            downloadAndMerge(ftp_link,dir)

        # if isComplete2020 and os.path.exists(out_file):
        #     shutil.copyfile(out_file, os.path.join(out_dir_complete,os.path.basename(out_file)))
print("Number completed samples up to 2020:"+str(num_complete_before_2020))
for y in range(2008,2023):
    if y in s.keys():
        print("Year "+str(y)+":"+str(s[y]["count"]))
        for m in range(1,13):
            if m in s[y].keys():
                print(" -Month "+str(m)+":"+str(s[y][m]))
threads=16
species="Escherichia coli"
def run_panta(input_dir, out_dir,log_file):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	cmd = ('/usr/bin/time -v python pan-genome.py main '
	  	   f'-o {out_dir} -t {threads} -g {input_dir}/*.gff '
		   f'1>> {log_file} 2>&1')
	os.system(cmd)
	return out_dir
def run_panta_add(basedir, gffdir,log_file):

	cmd = (f'/usr/bin/time -v python pan-genome.py main -o {basedir} -t {threads} -g {gffdir}/*.gff 1>> {log_file} 2>&1')
	os.system(cmd)
	return out_dir
def run_roary(input_dir, out_dir,log_file):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v roary -v -z -p {threads} -f {out_dir} {input_dir}/*.gff 1>> {log_file} 2>&1')
	os.system(cmd)

	#os.system(f'mv roary.log {log_file}')

	return out_dir
def run_pirate(input_dir, out_dir,logfile):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v PIRATE -i {input_dir} -o {out_dir} -t {threads} -z 2  1>> {logfile} 2>&1')
	os.system(cmd)

	#os.system(f'mv roary.log {log_file}')

	return out_dir

def run_panaroo(input_dir, out_dir,logfile):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v panaroo -i {input_dir}/*.gff -o {out_dir} -t {threads} --clean-mode strict --remove-invalid-genes 1>> {logfile} 2>&1')
	os.system(cmd)

	#os.system(f'mv panaroo.log {out_dir}')

	return out_dir
def run_panx(input_dir, species,logfile):
	cmd = (f'/usr/bin/time -v ./panX.py -fn {input_dir} -sl {species} '
		   f'-dmdc -sitr -t {threads} -cg 0.99 -slt -rlt 1>> {logfile} 2>&1')
	os.system(cmd)
	# os.system(f'mv panx.log {input_dir}')
	return input_dir
def backupfolder(sourcefolder,desfolder):
    cmd=f'cp -r {sourcefolder} {desfolder}'
    os.system(cmd)
#running section

log_folder=os.path.join(out_dir,"logs")
if not os.path.exists(log_folder):
    os.mkdir(log_folder)

# run complete genome folder:
incremental_gff=os.path.join(out_dir,"gffstore")
if not os.path.exists(incremental_gff):
    os.mkdir(incremental_gff)

cmd=f'cp -r {out_dir_complete}/*gff {incremental_gff}'
os.system(cmd)
run_dir_panta=os.path.join(out_dir,"panta_out")
run_panta(incremental_gff,run_dir_panta,os.path.join(log_folder,"pant_start.log"))
backupfolder(run_dir_panta,os.path.join(out_dir,"out_panta_start"))
run_dir_roary=os.path.join(out_dir,"roary_start")
run_roary(incremental_gff,run_dir_roary,os.path.join(log_folder,"roary_start.log"))
run_dir_pirate=os.path.join(out_dir,"pirate_start")
run_pirate(incremental_gff,run_dir_pirate,os.path.join(log_folder,"pirate_start.log"))
run_dir_panaroo=os.path.join(out_dir,"panaroo_start")
run_panaroo(incremental_gff,run_dir_panaroo,os.path.join(log_folder,"panaroo_start.log"))
#run_dir_panx=os.path.join(out_dir,"panx_start")
#run_panx(incremental_gff,run_dir_panx,os.path.join(log_folder,"panx_start.log"))
# run for each quater
for Q in ['Q1','Q2','Q3','Q4']:
    run_panta_add(run_dir_panta,os.path.join(out_dir,Q),os.path.join(log_folder,"pant_"+Q+".log"))
    backupfolder(run_dir_panta,os.path.join(out_dir,"out_panta_"+Q))
    #create incremental gff folder:
    cmd=f'cp -r {os.path.join(out_dir,Q)}/*gff {incremental_gff}'
    os.system(cmd)
    run_dir_roary=os.path.join(out_dir,"roary_"+Q)
    run_roary(incremental_gff,run_dir_roary,os.path.join(log_folder,"roary_"+Q+".log"))
    run_dir_pirate=os.path.join(out_dir,"pirate_"+Q)
    run_pirate(incremental_gff,run_dir_pirate,os.path.join(log_folder,"pirate_"+Q+".log"))
    run_dir_panaroo=os.path.join(out_dir,"panaroo_"+Q)
    run_panaroo(incremental_gff,run_dir_panaroo,os.path.join(log_folder,"panaroo_"+Q+".log"))
    #run_dir_panx=os.path.join(out_dir,"panx_"+Q)
    #run_panx(incremental_gff,run_dir_panx,os.path.join(log_folder,"panx_"+Q+".log"))
