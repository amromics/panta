import os
from glob import glob
import random
import shutil

def get_gff(source_dataset):
	gff_files = glob(source_dataset)
	gff_files.sort()
	random.seed(62)
	random.shuffle(gff_files)
	return gff_files


def copy_file(gff_list, out_dir):
	if not os.path.exists(out_dir):
		os.mkdir(out_dir)
	for gff in gff_list:
		os.system(f'cp {gff} {out_dir}')


def run_panta(input_dir, out_dir):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	
	cmd = ('/usr/bin/time -v python pan-genome.py main '
	  	   f'-o {out_dir} -t {threads} -g {input_dir}/*.gff '
		   f'1>> {out_dir}/panta.log 2>&1')
	os.system(cmd)

	return out_dir

def run_add_pipeline(gff_list, out_dir, n):
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	cmd = ('/usr/bin/time -v python pan-genome.py main '
		   f'-o {out_dir}/{str(n)} -t {threads} -g ')
	cmd += ' '.join(gff_list[:n])
	cmd += f' 1>> {out_dir}/{str(n)}.log 2>&1'
	os.system(cmd)
			
	cmd = ('/usr/bin/time -v python pan-genome.py add '
		   f'-o {out_dir}/{str(n)} -t {threads} -g ')
	cmd += ' '.join(gff_list[n:])
	cmd += f' 1>> {out_dir}/{str(n)}.log 2>&1'
	os.system(cmd)

	return out_dir

def run_roary(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)
	
	cmd = (f'/usr/bin/time -v roary -v -z -p {threads} -f {out_dir} '
		   f'{input_dir}/*.gff 1>> roary.log 2>&1')
	os.system(cmd)

	os.system(f'mv roary.log {out_dir}')

	return out_dir

def run_roary_nosplit(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v roary -v -s -z -p {threads} -f {out_dir} '
		   f'{input_dir}/*.gff 1>> roary.log 2>&1')
	os.system(cmd)

	os.system(f'mv roary.log {out_dir}')

	return out_dir


def run_roary_nosplit_align(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v roary -v -s -z -e -n -p {threads} -f {out_dir} '
		   f'{input_dir}/*.gff 1>> roary.log 2>&1')
	os.system(cmd)

	os.system(f'mv roary.log {out_dir}')

	return out_dir


def run_pirate(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v PIRATE -i {input_dir} -o {out_dir} '
		   f'-t {threads} -z 2  1>> pirate.log 2>&1')
	os.system(cmd)

	os.system(f'mv pirate.log {out_dir}')

	return out_dir


def run_pirate_align(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v PIRATE -i {input_dir} -o {out_dir} '
		   f'-t {threads} -z 2 --align  1>> pirate.log 2>&1')
	os.system(cmd)

	os.system(f'mv pirate.log {out_dir}')

	return out_dir

def run_panaroo(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v panaroo -i {input_dir}/*.gff -o {out_dir} '
		   f'-t {threads} --clean-mode strict 1>> panaroo.log 2>&1')
	os.system(cmd)

	os.system(f'mv panaroo.log {out_dir}')

	return out_dir

def run_panaroo_align(input_dir, out_dir):
	if os.path.exists(out_dir):
		shutil.rmtree(out_dir)

	cmd = (f'/usr/bin/time -v panaroo -i {input_dir}/*.gff -o {out_dir} '
		   f'-t {threads} --clean-mode strict -a pan --aligner mafft '
		   '1>> panaroo.log 2>&1')
	os.system(cmd)

	os.system(f'mv panaroo.log {out_dir}')

	return out_dir

def run_panx(input_dir, species):

	cmd = (f'/usr/bin/time -v ./panX.py -fn {input_dir} -sl {species} '
		   f'-dmdc -sitr -t {threads} -cg 0.99 -slt -rlt 1>> panx.log 2>&1')
	os.system(cmd)

	# os.system(f'mv panx.log {input_dir}')

	return input_dir



if __name__ == "__main__":
	base_dir = '/home/ntanh1999'
	collection_name = 'Pa100'
	panta_dir = f'{base_dir}/amromics/amromics/pan-genome'

	global threads
	threads = 8
	
	gff_dir = f'{base_dir}/{collection_name}/gff'
	
	# source_dataset = f'{base_dir}/data/Sp616_gff/*.gff'
	# gff_list = get_gff(source_dataset)
	# copy_file(gff_list[:100], gff_dir)
	

	# os.chdir(panta_dir)
	# run_panta(gff_dir, f'{base_dir}/{collection_name}/out/panta')
	# run_roary(gff_dir, f'{base_dir}/{collection_name}/out/roary')
	# run_roary_nosplit(gff_dir, f'{base_dir}/{collection_name}/out/roary_nosplit')
	# run_roary_nosplit_align(gff_dir, f'{base_dir}/{collection_name}/out/roary_nosplit_align')
	# run_pirate(gff_dir, f'{base_dir}/{collection_name}/out/PIRATE')
	# run_panaroo(gff_dir, f'{base_dir}/{collection_name}/out/panaroo')
	
	# os.chdir(f'{base_dir}/pan-genome-analysis')
	# run_panx(f'{base_dir}/{collection_name}/test', 'Klebsiella_pneumoniae')

	# run_pirate_align(gff_dir, f'{base_dir}/{collection_name}/out/PIRATE_align')
	# run_panaroo_align(gff_dir, f'{base_dir}/{collection_name}/out/panaroo_align')

	gff_list = get_gff(gff_dir + '/*.gff')
	os.chdir(panta_dir)
	for n in [2, 25, 50, 75, 99]:
		run_add_pipeline(
			gff_list, f'{base_dir}/{collection_name}/out/panta_add', n)