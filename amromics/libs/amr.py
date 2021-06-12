import os
import shutil
import csv
import logging
import multiprocessing
import gzip
import pandas as pd
import amromics.libs.bioseq as bioseq
import amromics.libs.element_finder as element_finder
from amromics.utils.command import run_command
import amromics.libs.mlst as mlst
logger = logging.getLogger(__name__)
NUM_CORES_DEFAULT = multiprocessing.cpu_count()
def detect_amr_abricate(prefix_name, assembly, base_dir='.', threads=8, overwrite=False, timing_log=None):
    """
    Run abricate to identify resistant genes

    Parameters
    ----------
    sample:
        a dictionary-like object containing various attributes for a sample
    sample_dir: str
        the directory of the sample
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to resistant gene file
    """
    path_out = os.path.join(base_dir,   prefix_name+'_abricate')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    # TODO: replace by consensus db later
    amr_out = os.path.join(path_out, prefix_name+ '_resistome.tsv')
    if os.path.isfile(amr_out) and (not overwrite):
        logger.info('Resistome for {} exists, skip analysis'.format(prefix_name))
        return amr_out
    dbs=['ncbi','megares','ecoh','argannot','card']
    numError=0
    outputfiles=[]
    for db in dbs:
        outfile= os.path.join(path_out,prefix_name + '_'+db+'.tsv')
        cmd = 'abricate --quiet --threads {threads} --nopath --db {db} {infile} > {outfile}'.format(
            threads=threads,
            db=db,
            infile=assembly,
            outfile=outfile)
        if run_command(cmd, timing_log) != 0:
            numError=numError+1
        else:
            outputfiles.append(outfile)
    if numError==len(dbs):
        raise Exception('Error running amr')
    combined_tsv = pd.concat([pd.read_csv(f,sep='\t') for f in outputfiles ])
    combined_tsv.sort_values(['SEQUENCE','START'],ascending=[True, True],inplace=True)
    combined_tsv.to_csv(amr_out, index=False,sep='\t', encoding='utf-8-sig')
    #sample['updated'] = True
    return amr_out
def detect_amr(prefix_name,faa_file,fna_file,gff_file,genus=None,species=None,  base_dir='.', db='db/amrfinderplus/data/latest', timing_log=None, threads=0):
    """
        Run AMR analysis, using AMRfinderPlus for searching amr genes, virulome genes and point mutaions.
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT
    path_out = os.path.join(base_dir,   prefix_name+'_amrfinder')
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    #AMR profiling with CARD. TODO: replace by consensus db later
    ret_out = os.path.join(path_out, prefix_name + '_amr.tsv')
    amr_out = os.path.join(path_out, prefix_name+ '_resistome.tsv')
    virulen_out = os.path.join(path_out, prefix_name + '_virulome.tsv')
    point_out = os.path.join(path_out, prefix_name + '_point.tsv')
    #using abricate
    # cmd = 'abricate --quiet --threads {threads} --nopath --db card {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=amr_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    #using build-in function
    #element_finder.search_amr(sample=read_data['assembly'],output=amr_out,threads=threads)
    #using AMRFinderPlus

    #process files in prokka folder, prepare for using amrfinder
    #move files from prokka to temp folder
    temp_dir = os.path.join(base_dir, 'amr_temp')
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir)

    temp_gff_file=os.path.join(temp_dir,prefix_name+'.gff')
    source_gff_file=None
    # for root, dirs, files in os.walk(read_data['annotation']):
    #     for _file in files:
    #         if _file.endswith(('.faa')):
    #             faa_file = shutil.copyfile(os.path.join(str(root),_file), faa_file)
    #         if _file.endswith(('.fna')):
    #             fna_file = shutil.copyfile(os.path.join(str(root),_file), fna_file)
    #         if _file.endswith(('.gff')):
    #             source_gff_file = os.path.join(str(root),_file)

    source_gff_file=gff_file
    #add Name property to column 9 of gff file (AMRfinder need it!) and remove #fasta section
    if not source_gff_file==None:
        destination= open(temp_gff_file, "w" )
        #source= open( source_gff_file, "r" )
        with gzip.open(source_gff_file,'rt') as source:
            for line in source:
                if line.startswith('##FASTA'):
                    break
                if line.startswith('##'):
                    destination.write( line )
                else:
                    newline=line.replace('ID=','Name=')
                    destination.write( newline )
        #source.close()
        destination.close()
    gunzip_faa= faa_file;
    if faa_file.endswith('.gz'):
        gunzip_faa =os.path.join(temp_dir,prefix_name+'.faa')
        cmd = 'gunzip -c {} > {}'.format(faa_file, gunzip_faa)
        run_command(cmd)
    gunzip_fna= fna_file;
    if fna_file.endswith('.gz'):
        gunzip_fna =os.path.join(temp_dir,prefix_name+'.fna')
        cmd = 'gunzip -c {} > {}'.format(fna_file, gunzip_fna)
        run_command(cmd)
    cmd = 'amrfinder -d {database} -p {faa_file}  -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile}'\
    .format(
        database=db,
        faa_file=gunzip_faa,
        fna_file=gunzip_fna,
        gff_file=temp_gff_file,
        threads=threads,
        outfile=ret_out
    )
    #full option if has --Genus
    if not genus==None:
        organism = genus.capitalize()
        if not species==None and not species=='':
            organism = species.replace(' ','_')
        organisms = ['Campylobacter', 'Enterococcus_faecalis', 'Enterococcus_faecium', 'Escherichia', 'Klebsiella', 'Salmonella', 'Staphylococcus_aureus', 'Staphylococcus_pseudintermedius', 'Vibrio_cholerae']
        if organism in organisms:
            cmd = 'amrfinder -d {database} -p {faa_file} -O {organism}  -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile}'\
            .format(
                database=db,
                faa_file=faa_file,
                organism=organism,
                fna_file=fna_file,
                gff_file=gff_file,
                threads=threads,
                outfile=ret_out
            )
        elif read_data['genus'] in organisms:
            cmd = 'amrfinder -d {database} -p {faa_file} -O {genus}  -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile}'\
            .format(
                database=db,
                faa_file=faa_file,
                genus=read_data['genus'],
                fna_file=fna_file,
                gff_file=gff_file,
                threads=threads,
                outfile=ret_out
            )
        else:
            cmd = 'amrfinder -d {database} -p {faa_file} -n {fna_file} -g {gff_file} --plus --threads {threads} -o {outfile}'\
            .format(
                database=db,
                faa_file=faa_file,
                fna_file=fna_file,
                gff_file=gff_file,
                threads=threads,
                outfile=ret_out
            )

    cmd = "bash -c '{}'".format(cmd)
    if run_command(cmd,timing_log) != 0:
        return None

    #clean up:
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    #proccess output files:
    virulen=[]
    amr=[]
    point=[]
    header=[]
    with open(ret_out) as tsvfile:
        reader = csv.DictReader(tsvfile, dialect='excel-tab')
        for row in reader:
            header=row.keys()
            if row['Element type']=='VIRULENCE':
                virulen.append(row)
            elif row['Element subtype']=='POINT':
                point.append(row)
            else:
                amr.append(row)
    with open(amr_out, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=header,delimiter='\t')
        writer.writeheader()
        for row in amr:
            writer.writerow(row)

    with open(point_out, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=header,delimiter='\t')
        writer.writeheader()
        for row in point:
            writer.writerow(row)

    with open(virulen_out, 'w', newline='') as csvfile:

        writer = csv.DictWriter(csvfile, fieldnames=header,delimiter='\t')
        writer.writeheader()
        for row in virulen:
            writer.writerow(row)

    if os.path.exists(ret_out):
        os.remove(ret_out)

    return amr_out,point_out,virulen_out

###Virulome profiling using abricate with VFDB
def detect_virulome(prefix_name,assembly, base_dir='.', threads=0, timing_log=None):
    """
    Run in-house script to identify virulent genes using VFDB

    Parameters
    ----------
    prefix_name:
        name to attach to output
    assembly: str
        input sequence
    threads: int
        number of threads to use
    overwrite:bool
        whether to overwrite the existing result
    timing_log: str
        log file
    Returns
    -------
        path to virulent gene file
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir,  prefix_name+'_element_finder' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    vir_out = os.path.join(path_out, prefix_name + '_virulome.tsv')
    # cmd = 'abricate --quiet --threads {threads} --nopath --db vfdb {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=vir_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    gunzip_fna= assembly;
    if assembly.endswith('.gz'):
        gunzip_fna =os.path.join(path_out,prefix_name+'.fasta')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fna)
        run_command(cmd)
    element_finder.search_virulome(sample=gunzip_fna,output=vir_out,threads=threads)

    if not os.path.exists(os.path.join(path_out,prefix_name+'.fasta')):
        os.remove(os.path.join(path_out,prefix_name+'.fasta'))
    return vir_out

###Find plasmid's origin of replication using abricate with plasmidfinder db
def detect_plasmid(prefix_name,assembly,  base_dir='.', threads=0, timing_log=None):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir,  prefix_name+'_element_finder')
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    #Plasmid finder
    oriREP_out = os.path.join(path_out,prefix_name + '_plasmid.tsv')
    # cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=oriREP_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    gunzip_fna= assembly;
    if assembly.endswith('.gz'):
        gunzip_fna =os.path.join(path_out,prefix_name+'.fasta')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fna)
        run_command(cmd)
    element_finder.search_plasmid(sample=assembly,output=oriREP_out,threads=threads)
    if not os.path.exists(os.path.join(path_out,prefix_name+'.fasta')):
        os.remove(os.path.join(path_out,prefix_name+'.fasta'))

    return oriREP_out
def detect_pmlst(prefix_name,assembly,  base_dir='.', threads=0):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir,prefix_name+'_pmlst' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    #Plasmid finder
    pmlst_out = os.path.join(path_out, prefix_name + '_pmlst.tsv')
    # cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=oriREP_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    gunzip_fna= assembly;
    if assembly.endswith('.gz'):
        gunzip_fna =os.path.join(path_out,prefix_name+'.fasta')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fna)
        run_command(cmd)
    m=mlst.find_mlst(query_file=gunzip_fna,blastdb='db/pmlst/blast/pmlst.fa',mlstdb='db/pmlst/pubmlst',num_threads=threads)
    with open(pmlst_out, 'w') as f:
        f.write("%s\t%s\t%s"%(m['file'],m['scheme'],m['st']))
        for gene in m['profile']:
            f.write("\t%s"%gene)
        f.write("\n")
    # cmd = 'mlst --quiet --threads {threads} --nopath {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=mlst_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    if os.path.exists(os.path.join(path_out,prefix_name+'.fasta')):
        os.remove(os.path.join(path_out,prefix_name+'.fasta'))



    return pmlst_out
def detect_insertion_sequence(prefix_name,assembly,  base_dir='.', threads=0):
    """
        Run isescan for searching IS
        :param read_data: result holder
        :param base_dir: working directory
        :return: path to output file in result holder
    """
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, prefix_name+"_isescan")
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    gunzip_fna= assembly;
    if assembly.endswith('.gz'):
        gunzip_fna =os.path.join(path_out,prefix_name+'.fasta')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fna)
        run_command(cmd)
    #Plasmid finder
    isescan_out = os.path.join(path_out, prefix_name + '_is.tsv')
    cmd = 'isescan.py --nthread {threads} --seqfile {asm} --output {output}  '.format(
        threads=threads,

        asm=gunzip_fna,
        output=path_out


    )
    if run_command(cmd) != 0:
        return None
    #if os.path.exists(path_out+'/prediction'):
    #    shutil.rmtree(path_out+'/prediction')
    #shutil.copytree('prediction', path_out+'/prediction')
    #read_data['is'] = path_out+'/prediction'
    isout=None
    for root, dirs, files in os.walk(path_out, topdown=False):
        for name in files:
            if name.endswith('.raw'):
                isout=os.path.join(root, name)
    #if os.path.exists('prediction'):
    #    shutil.rmtree('prediction')
    if os.path.exists(os.path.join(path_out,prefix_name+'.fasta')):
        os.remove(os.path.join(path_out,prefix_name+'.fasta'))
    return isout

def detect_integron(prefix_name,assembly, base_dir='.', timing_log=None,threads=0):
    if threads == 0:
        threads = NUM_CORES_DEFAULT
    #path_out = os.path.join(base_dir, 'integron_finder_' + read_data['sample_id'])
    path_out = os.path.join(base_dir,  prefix_name+'_integrall' )
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    #Plasmid finder
    gunzip_fna= assembly;
    if assembly.endswith('.gz'):
        gunzip_fna =os.path.join(path_out,prefix_name+'.fasta')
        cmd = 'gunzip -c {} > {}'.format(assembly, gunzip_fna)
        run_command(cmd)
    integron_out = os.path.join(path_out,prefix_name + '_integron.tsv')
    # cmd = 'integron_finder {sequence} --func-annot --local-max --mute --outdir {outdir}'.format(sequence=read_data['assembly'],outdir=path_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd,timing_log) != 0:
    #     return None
    element_finder.search_integrall(sample=gunzip_fna,output=integron_out,threads=threads)
    if os.path.exists(os.path.join(path_out,prefix_name+'.fasta')):
        os.remove(os.path.join(path_out,prefix_name+'.fasta'))
    return integron_out
def detect_prophage(prefix_name,faa_file, base_dir='.', timing_log=None,threads=0):
    if threads == 0:
        threads = NUM_CORES_DEFAULT

    path_out = os.path.join(base_dir, 'element_finder_' +prefix_name)
    if not os.path.exists(path_out):
        os.makedirs(path_out)

    #Plasmid finder
    prophage_out = os.path.join(path_out,prefix_name + '_prophage.tsv')
    # cmd = 'abricate --quiet --threads {threads} --nopath --db plasmidfinder {infile} > {outfile}'.format(threads=threads,infile=read_data['assembly'],outfile=oriREP_out)
    # cmd = "bash -c '{}'".format(cmd)
    # if run_command(cmd) != 0:
    #     return None
    #Plasmid finder
    gunzip_faa= faa_file;
    if faa_file.endswith('.gz'):
        gunzip_faa =os.path.join(path_out,prefix_name+'.faa')
        cmd = 'gunzip -c {} > {}'.format(faa_file, gunzip_faa)
        run_command(cmd)
    element_finder.search_prophage(sample=gunzip_faa,output=prophage_out,threads=threads)
    if os.path.exists(os.path.join(path_out,prefix_name+'.faa')):
        os.remove(os.path.join(path_out,prefix_name+'.faa'))

    return prophage_out
