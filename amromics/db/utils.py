# -*- coding: utf-8 -*-
"""
    Utilities


Revision history:
----------------
2019-08-16: MDC Created

"""
from __future__ import division, print_function, absolute_import

import logging
import os, glob,shutil
import amromics.libs.bioseq as bioseq
from collections import OrderedDict
import xml.etree.ElementTree as ET
import wget
import urllib
import gzip
import shutil
from zipfile import ZipFile
import re

logger = logging.getLogger(__name__)



def downloadfile(url,savefile):
    #save the bak file before download new file if file exists
    if os.path.exists(savefile):
        os.rename(savefile,savefile+'.bak')
    try:
        print('\nDownloading file {}'.format(savefile))
        wget.download(url,savefile)
        if os.path.exists(savefile+'.bak'):
            os.remove(savefile+'.bak')
        return savefile
    except BaseException as error:
        print(str(error))
        if os.path.exists(savefile+'.bak'):
            os.rename(savefile+'.bak',savefile)
        return None
def get_mlst_db():
    '''
        Download PubMLST from ftp folder, read dbases.xml to get url for downloading
        combine sequences into one .fa file and make blast db
    '''
    url = 'https://pubmlst.org/data/dbases.xml'
    webdata = urllib.request.urlopen(url)
    tree = ET.parse(webdata)
    root = tree.getroot()
    mlst_dir='db/mlst'
    pubmlst_dir=os.path.join(mlst_dir,'pubmlst')
    if not os.path.exists(pubmlst_dir):
        os.makedirs(pubmlst_dir)

    for species in root:
        scheme_dir=''
        for profiles in species.iter('profiles'):
            for child in profiles:
                if child.tag=='url':
                    #scheme=os.path.basename(child.text).strip().replace('.txt','')
                    toks = child.text.split("/")
                    scheme_name = toks[toks.index("db")+1].split("_")[1]
                    scheme_num = toks[toks.index("schemes")+1]
                    scheme = scheme_name + (("_"+scheme_num) if scheme_num!="1" else "")
                    scheme_dir=os.path.join(pubmlst_dir,scheme)
                    if not os.path.exists(scheme_dir):
                        os.makedirs(scheme_dir)
                    #print(os.path.basename(child.text).strip())
                    downloadfile(child.text,os.path.join(scheme_dir,scheme+".txt"))
        for loci in species.iter('loci'):
            for locus in loci:
                for child in locus:
                    if child.tag=='url':
                        print(os.path.join(scheme_dir,locus.text.strip()))
                        downloadfile(child.text,os.path.join(scheme_dir,locus.text.strip()+".tfa"))
    ignore_species={'afumigatus','blastocystis','calbicans','cglabrata','ckrusei','ctropicalis','csinensis','kseptempunctata','sparasitica','tvaginalis'}
    for rm_species in ignore_species:
        rm_dir=os.path.join(pubmlst_dir,rm_species)
        if os.path.exists(rm_dir):
            shutil.rmtree(rm_dir)
    #make blastdb:
    blast_dir=os.path.join(mlst_dir,'blast')
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    mlst_file=os.path.join(blast_dir,'mlst.fa')
    try:
        os.remove(mlst_file)
    except OSError:
        pass
    for root, dirs, files in os.walk(pubmlst_dir):
        for _file in files:
            if _file.endswith('.tfa'):
                infile=str(root)+'/'+_file
                scheme = os.path.basename(str(root))
                catcmd='cat {tfa} | sed -e \'s/>/>{scheme}./g\' >>{mlst}'.format(
                    scheme=scheme,
                    tfa=infile,
                    mlst=mlst_file

                )
                #print(catcmd)
                os.system(catcmd)
    cmd="makeblastdb -hash_index -in {blastfile} -dbtype nucl -title PubMLST -parse_seqids".format(
        blastfile=mlst_file
    )
    print (cmd)
    os.system(cmd)
def get_virulome_db():
    '''
        Download vfdb from mgc.ac.cn, rewrite sequence'name
    '''
    vfdb_dir='db/vfdb'
    if not os.path.exists(vfdb_dir):
        os.makedirs(vfdb_dir)
    vfdb_file_zip=os.path.join(vfdb_dir,'sequences.gz')
    vfdb_file=os.path.join(vfdb_dir,'sequences.fa')
    vfdb_file_out=os.path.join(vfdb_dir,'sequences')
    vfdb_file_zip=downloadfile('http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz',vfdb_file_zip)
    if not vfdb_file_zip==None:
        with gzip.open(vfdb_file_zip, 'rb') as f_in:
            with open(vfdb_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        with open(vfdb_file_out,'w') as f:
            for seq in bioseq.read_sequence_file(vfdb_file):
                #print(seq.name+';'+seq.desc)
                accession=''
                z=re.match(r"^(\w+)\(\w+\|(\w+)(\.\d+)?\)", seq.name)
                if z:
                    accession=z[2]
                else:
                    accession=seq.name
                z=re.match(r"^\((.*?)\)", seq.desc)
                gene=''
                if z:
                    gene=z[1]
                else:
                    gene=accession
                seq.set_name('vfdb~~~'+gene+'~~~'+accession)
                f.write(seq.format_fasta(max_line=100))
    #clean up
    if os.path.exists(vfdb_file_zip):
        os.remove(vfdb_file_zip)
    if os.path.exists(vfdb_file):
        os.remove(vfdb_file)
    #make blast db
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title vfdb".format(
        blastfile=vfdb_file_out
    )
    print (cmd)
    os.system(cmd)
def get_plasmidfinder():
    plasmidfinder_dir='db/plasmidfinder'
    if not os.path.exists(plasmidfinder_dir):
        os.makedirs(plasmidfinder_dir)
    url='https://bitbucket.org/genomicepidemiology/plasmidfinder_db/get/HEAD.zip'
    plasmidfinder_zip=os.path.join(plasmidfinder_dir,'plasmidfinder.zip')
    plasmidfinder_unzip=os.path.join(plasmidfinder_dir,'temp')
    plasmidfinder_out=os.path.join(plasmidfinder_dir,'sequences')
    plasmidfinder_zip=downloadfile(url,plasmidfinder_zip)
    if not plasmidfinder_zip==None:

        with ZipFile(plasmidfinder_zip, 'r') as zipObj:
            zipObj.extractall(plasmidfinder_unzip)
        with open(plasmidfinder_out,'w') as f:
            for root, dirs, files in os.walk(plasmidfinder_unzip):
                for _file in files:
                    if _file.endswith(('.fsa')):
                        #print(str(root)+'/'+_file)
                        for seq in bioseq.read_sequence_file(str(root)+'/'+_file):
                            accession=''
                            gene=''
                            #parse string to get gene and acc
                            z=re.match(r"^(.*)_(([A-Z]+|NC_)\d+(\.\d+)?)", seq.name)
                            if z:
                                accession=z[2]
                                gene=z[1]
                            else:
                                accession=seq.name

                            seq.set_desc(seq.name)
                            seq.set_name('plasmidfinder~~~'+gene+'~~~'+accession+'~~~')
                            f.write(seq.format_fasta(max_line=100))
    #cleanup
    if os.path.exists(plasmidfinder_zip):
        os.remove(plasmidfinder_zip)
    if os.path.exists(plasmidfinder_unzip):
        shutil.rmtree(plasmidfinder_unzip)
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title plasmidfinder".format(
        blastfile=plasmidfinder_out
    )
    print (cmd)
    os.system(cmd)
def get_pmlst():
    '''
    Get db from ge.cbs.dtu.dk, format db like https://github.com/tseemann/mlst
    '''
    pmlst_dir='db/pmlst'
    if not os.path.exists(pmlst_dir):
        os.makedirs(pmlst_dir)
    blast_dir=os.path.join(pmlst_dir,'blast')
    pubmlst_dir=os.path.join(pmlst_dir,'pubmlst')
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    url='https://bitbucket.org/genomicepidemiology/pmlst_db/get/HEAD.zip'
    pmlst_zip=os.path.join(pmlst_dir,'pmlst.zip')
    pmlst_unzip=os.path.join(pmlst_dir,'temp')

    pmlst_zip=downloadfile(url,pmlst_zip)
    pmlst_file=os.path.join(blast_dir,'pmlst.fa')
    try:
        os.remove(pmlst_file)
    except OSError:
        pass
    if not pmlst_zip==None:
        with ZipFile(pmlst_zip, 'r') as zipObj:
            zipObj.extractall(pmlst_unzip)
        for root, dirs, files in os.walk(pmlst_unzip):
            for _file in files:
                if _file.endswith(('.fsa')):
                    #print(str(root)+'/'+_file)
                    infile=str(root)+'/'+_file
                    scheme = os.path.basename(_file).strip().replace('.fsa','')
                    catcmd='cat {fsa} | sed -e \'s/>/>{scheme}./g\' >>{pmlst}'.format(
                        scheme=scheme,
                        fsa=infile,
                        pmlst=pmlst_file

                    )
                    #print(catcmd)
                    os.system(catcmd)
                if _file.endswith(('.txt.clean')):
                    scheme=os.path.basename(_file).strip().replace('.txt.clean','')
                    scheme_dir=os.path.join(pubmlst_dir,scheme)
                    if not os.path.exists(scheme_dir):
                        os.makedirs(scheme_dir)
                    #move profile .txt to scheme folder
                    shutil.move(os.path.join(root,_file), os.path.join(scheme_dir,scheme+'.txt'))
    #cleanup
    if os.path.exists(pmlst_zip):
        os.remove(pmlst_zip)
    if os.path.exists(pmlst_unzip):
        shutil.rmtree(pmlst_unzip)
    cmd="makeblastdb -hash_index -in {blastfile} -dbtype nucl -title pMLST -parse_seqids".format(
        blastfile=pmlst_file
    )
    print (cmd)
    os.system(cmd)
    #https://bitbucket.org/genomicepidemiology/intfinder_db
def get_integron():
    integron_dir='db/integron'
    if not os.path.exists(integron_dir):
        os.makedirs(integron_dir)
    url='https://bitbucket.org/genomicepidemiology/intfinder_db/get/HEAD.zip'
    integron_zip=os.path.join(integron_dir,'integron.zip')
    integron_unzip=os.path.join(integron_dir,'temp')
    integron_out=os.path.join(integron_dir,'sequences')
    integron_zip=downloadfile(url,integron_zip)
    if not integron_zip==None:

        with ZipFile(integron_zip, 'r') as zipObj:
            zipObj.extractall(integron_unzip)
        with open(integron_out,'w') as f:
            for root, dirs, files in os.walk(integron_unzip):
                for _file in files:
                    if _file.endswith(('.fsa')):
                        #print(str(root)+'/'+_file)
                        for seq in bioseq.read_sequence_file(str(root)+'/'+_file):
                            accession=''
                            gene=''
                            #parse string to get gene and acc
                            z=re.match(r"^(.*)_(([A-Z]+|NC_)\d+(\.\d+)?)", seq.name)
                            if z:
                                accession=z[2]
                                gene=z[1]
                            else:
                                accession=seq.name

                            seq.set_desc(seq.name)
                            seq.set_name('integron~~~'+gene+'~~~'+accession+'~~~')
                            f.write(seq.format_fasta(max_line=100))
    #cleanup
    if os.path.exists(integron_zip):
        os.remove(integron_zip)
    if os.path.exists(integron_unzip):
        shutil.rmtree(integron_unzip)
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title integron".format(
        blastfile=integron_out
    )
    print (cmd)
    os.system(cmd)
def get_kraken2():
    '''
        Get kraken db
    '''
    #url='https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20200919.tar.gz'
    #for better download speed
    url='ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/old/minikraken2_v1_8GB_201904.tgz'
    kraken2_db='db/kraken2/'
    if not os.path.exists(kraken2_db):
        os.makedirs(kraken2_db)
    k2std_zip=os.path.join(kraken2_db,'k2std.tar.gz')
    k2std_unzip=os.path.join(kraken2_db,'k2std')
    k2std_zip=downloadfile(url,k2std_zip)
    if not k2std_zip==None:
        shutil.unpack_archive(k2std_zip, k2std_unzip)
        for (root,dirs,files) in os.walk(k2std_unzip, topdown=True):
            for f in files:
                filepath=os.path.join(root, f)
                shutil.move(filepath, os.path.join(k2std_unzip,f))
    if os.path.exists(k2std_zip):
        os.remove(k2std_zip)
def get_prophage():
    #url='http://phaster.ca/downloads/prophage_virus.db'
    url='https://gdurl.com/y2eu'
    prophage_dir='db/prophage'
    if not os.path.exists(prophage_dir):
        os.makedirs(prophage_dir)
    prophage_file=os.path.join(prophage_dir,'sequences')
    prophage_file=downloadfile(url,prophage_file)
    with open(integron_out,'w') as f:
        for seq in bioseq.read_sequence_file(prophage_file):
            accession=''
            gene=''
            #parse string to get gene and acc
            z=re.match(r"^PROPHAGE_(.*)\|(.*)\|(.*)\|(([A-Z]+|NP_)\d+(\.\d+)?)", seq.name)
            if z:
                accession=z[4]
                gene=z[1]
            else:
                accession=seq.name

                seq.set_desc(seq.name)
            seq.set_name('integron~~~'+gene+'~~~'+accession+'~~~')
            f.write(seq.format_fasta(max_line=100))
    cmd="makeblastdb -in {blastfile} -dbtype prot -title prophage".format(
        blastfile=prophage_file
    )
    print (cmd)
    os.system(cmd)
def update_amrfinderplus():
    '''
        Get amrfinderplus db and save to custom folder
    '''
    amrfinderplus_dir='db/amrfinderplus/data/'
    if not os.path.exists(amrfinderplus_dir):
        os.makedirs(amrfinderplus_dir)

    cmd="amrfinder_update -d {amrfinderplusdb}".format(amrfinderplusdb=amrfinderplus_dir)
    os.system(cmd)

def get_integrall_db():
    integrall_dir='db/integrall'
    if not os.path.exists(integrall_dir):
        os.makedirs(integrall_dir)
    url='http://integrall.bio.ua.pt/?getFastaAll'
    integrall_file=os.path.join(integrall_dir,'integrall.fasta')
    integrall_file_out=os.path.join(integrall_dir,'sequences')
    integrall_file=downloadfile(url,integrall_file)
    if not integrall_file==None:
        with open(integrall_file_out,'w') as f:
            for seq in bioseq.read_sequence_file(integrall_file):
                print(seq.name)
                if '?' in seq.name:
                    continue
                z= seq.name.split('|')
                y=seq.desc.split('|')
                if 'Δ' in y[1]:
                    y[1]=y[1].replace('Δ','beta-')
                seq.set_name('integrall~~~'+y[1]+'~~~'+z[1])
                seq.set_desc(y[0])
                f.write(seq.format_fasta(max_line=100))
    cmd="makeblastdb -in {blastfile} -dbtype nucl -title integrall".format(
        blastfile=integrall_file_out
    )
    print (cmd)
    os.system(cmd)
def setup_db():
    get_kraken2()
    get_mlst_db()
    get_virulome_db()
    get_plasmidfinder()
    update_amrfinderplus()
    get_pmlst()
    get_integrall_db()
    get_integron()
    #get_prophage()
def setup_minidb():
    get_mlst_db()
    get_virulome_db()
