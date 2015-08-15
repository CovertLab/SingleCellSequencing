"""
Import python packages
"""

import HTSeq 
import collections
import itertools
import os
import subprocess
import collections
import datetime
import yaml
import fnmatch
import shlex
import numpy
import scipy
import scipy.io as sio 

"""
Define functions
"""

def run_cmd(cmd):
	environ = {
		"PATH": os.environ["PATH"],
		"LANG": "C",
		"LC_ALL": "C",
	}
	out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env = environ).communicate()[0]
	return out

def write_file(filename, content):
	h = open(filename, "w")
	h.write(content)
	h.close()

def parse_filename(filename):
	cell_ID = filename.split('_')[0]
	return cell_ID

def make_directories(direc_name, subdirecs_to_make = None):
	for direc in subdirecs_to_make:
		if not os.path.exists(os.path.join(direc_name, direc)):
			os.makedirs(os.path.join(direc_name, direc))

def unzip_file(filename, input_direc = None, output_direc = None):
	cmd = ['gunzip', input_direc + filename]
	run_cmd(cmd)
	cmd = ['cp', input_direc + filename[:-3], output_direc + filename[:-3]]
	run_cmd(cmd)

def unzip_rawdata(direc_name):
	raw_data_direc = direc_name + 'Raw_Data/'
	unzip_direc = direc_name + 'unzipped/'
	file_list = os.listdir(raw_data_direc)
	for seq_file in file_list:
		print seq_file
		if fnmatch.fnmatch(seq_file, r'*.gz'):
			if not os.path.exists(unzip_direc + seq_file[:,-3]):
				unzip_file(seq_file, input_direc = raw_data_direc, output_direc = unzip_direc)

def trim_reads(direc_name, input_filename, unzipped_name = 'unzipped', trimmed_name = 'trimmed' ):
	
	left_file = os.path.join(direc_name, unzipped_name, input_filename + '_R1_001.fastq')
	right_file = os.path.join(direc_name, unzipped_name, input_filename + '_R2_001.fastq')


	ofn_temp = input_filename.split('_')
	output_filename = ofn_temp[0] + '_' + ofn_temp[1]

	output_left = os.path.join(direc_name, trimmed_name, output_filename + '_R1_trimmed.fq')
	output_right = os.path.join(direc_name, trimmed_name, output_filename + '_R2_trimmed.fq')

	if not os.path.exists(output_left):
		if not os.path.exists(output_right):
			nextera_file = '/share/PI/mcovert/downloads/bbmap_34.33/resources/nextera.fa.gz'
			trueseq_file = '/share/PI/mcovert/downloads/bbmap_34.33/resources/truseq.fa.gz'
			bbduk_params = ' overwrite=true ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=rl trimq=10 maq=15 minlen=50 ftm=5'

			cmd = 'bbduk.sh' + ' -Xmx2g' + ' in=' + left_file + ' in2=' + right_file + ' ref=' + nextera_file + ',' + trueseq_file + bbduk_params + ' out=' + output_left + ' out2=' + output_right

			cmd_temp = shlex.split(cmd)
			run_cmd(cmd_temp)

def trim_direc(direc_name):
	unzip_direc = direc_name + 'unzipped/'
	trimmed_direc = direc_name + 'trimmed/'
	file_list = os.listdir(unzip_direc)
	for seq_file in file_list:
		sft = seq_file.split('.')[0]
		sft2 = sft.split('_')
		seq_file_name = sft2[0] + '_' + sft2[1] + '_' + sft2[2] 
		print seq_file_name
		trim_reads(direc_name, seq_file_name)

def generate_genome_index():
	option_name = []
	option = []

	option_name += option_name + ['--runThreadN']
	option += ['16']
	option_name += ['--runMode'] 
	option += ['genomeGenerate']
	option_name += ['--genomeDir']
	option += ['/scratch/users/dvanva/Sequencing/Ref_seqs/']
	option_name += ['--genomeFastaFiles']
	option += ['/scratch/users/dvanva/Sequencing/Ref_seqs/mouse_withSpikeIns.fa']
	option_name += ['--sjdbGTFfile']
	option += ['/scratch/users/dvanva/Sequencing/Ref_seqs/mouse_withSpikeIns.gtf']
	option_name += ['--sjdbOverhang ']
	option += ['76']

	cmd = ['STAR']
	for j in xrange(len(option_name)):
		cmd = cmd + [option_name[j]] + [option[j]]

	print cmd
	run_cmd(cmd)

def run_STAR(filename):
	option_name = []
	option = []

	option_name[0] = ' --runThreadN '
	option[0] = '8'

def run_kallisto_index(input_filename, output_filename):
	cmd = ['kallisto', 'index', '-i', output_filename, input_filename]
	run_cmd(cmd)

def kallisto_direc(direc_name):
	output_directory = os.path.join(direc_name, 'align')
	file_list = os.listdir(os.path.join(direc_name,'trimmed'))
	for seq_file in file_list:
		sft = seq_file.split('.')[0]
		sft2 = sft.split('_')
		seq_file_name = sft[0] + '_' + sft2[1]
		if not os.path.exists(os.path.join(output_directory, seq_file_name + '.sam')):
			run_kallisto_quant(direc_name, seq_file_name)

def run_kallisto_quant(direc_name, input_filename, index_filename = '/scratch/PI/mcovert/dvanva/sequencing/reference_sequences/mus_musculus_ensembl_r81_index', trimmed_name = 'trimmed', aligned_name = 'aligned'): #, output_directory, output_filename, pair1, pair2):
	input_directory = os.path.join(direc_name,trimmed_name)
	output_directory = os.path.join(direc_name, aligned_name)

	pair1 = os.path.join(input_directory, input_filename+'_R1_trimmed.fq')
	pair2 = os.path.join(input_directory, input_filename+'_R2_trimmed.fq')
	cmd = ['kallisto', 'quant', '-i', index_filename, '-o', output_directory, '--pseudobam', pair1, pair2]
	
	# Use kallisto to align reads to the transcriptome and output the pseudoalignment to a sam file
	kal_out = run_cmd(cmd)
	write_file(os.path.join(output_directory, input_filename + '.sam'), kal_out)
	
	# Move files to the output directory
	h5_name = os.path.join(output_directory, input_filename+'.h5')
	cmd = ['mv', os.path.join(output_directory, 'abundance.h5'), h5_name]
	run_cmd(cmd)
	
	tsv_name = os.path.join(output_directory, input_filename+'.tsv')
	cmd = ['mv', os.path.join(output_directory, 'abundance.tsv'), tsv_name]
	run_cmd(cmd)

	run_info_name = os.path.join(output_directory, input_filename+'_runinfo.json')
	cmd = ['mv', os.path.join(output_directory, 'run_info.json'), run_info_name]
	run_cmd(cmd)

def run_kallisto_test():
	direc_name = '/scratch/PI/mcovert/dvanva/sequencing/library1'
	input_filename = '1-10_S14'
	run_kallisto_quant(direc_name, input_filename)

	# cmd = 'kallisto quant -i /scratch/PI/mcovert/dvanva/sequencing/reference_sequences/mus_musculus_index -o /scratch/PI/mcovert/dvanva/sequencing/library1 --pseudobam /scratch/PI/mcovert/dvanva/sequencing/library1/trimmed/1-10_S14_R1_trimmed.fq /scratch/PI/mcovert/dvanva/sequencing/library1/trimmed/1-10_S14_R2_trimmed.fq'
	# cmd_temp = shlex.split(cmd)
	# print cmd_temp
	# kal_out = run_cmd(cmd_temp)
	# write_file('/scratch/PI/mcovert/dvanva/sequencing/kal_out.sam', kal_out)

def load_sequence_counts(samfile_name = None, mouse_gtf = None, spikein_gtf = None):

	"""
	Load gtf files
	"""

	mouse_gtf_file = HTSeq.GFF_Reader(mouse_gtf)
	spikein_gtf_file = HTSeq.GFF_Reader(spikein_gtf)

	"""
	Load SAM file
	"""
	samfile = HTSeq.BAM_Reader(samfile_name)

	"""
	Construct list of genes for the mouse genome and spikeins
	"""
	
	spikein_transcript_list = []
	for feature in spikein_gtf_file:
		print feature
		if feature.type == 'exon':
			spikein_transcript_list += [feature.attr['gene_id']]

	print spikein_transcript_list

	mouse_exon_list = []
	for feature in mouse_gtf_file:
		if feature.type == 'exon':
			mouse_exon_list += [feature.attr['gene_id']]


	"""
	Construct genomic array for exons and spikeins
	"""
	exons = HTSeq.GenomicArrayOfSets('auto', stranded = False)
	spikeins = HTSeq.GenomicArrayOfSets('auto', stranded = False)

	for feature in mouse_gtf_file:
		if feature.type == 'exon':
			exons[ feature.iv ] += feature.attr['gene_id']

	for feature in spikein_gtf_file:
		if feature.type == 'exon':
			spikeins[ feature.iv ] += feature.attr['gene_id']

	"""
	Compute coverage map
	"""
	coverage = HTSeq.GenomicArray('auto', stranded = False)
	for bundle in HTSeq.pair_SAM_alignments(samfile, bundle = True):
		# Skip multiple alignments
		if len(bundle) != 1:
			continue 
		left_read, right_read = bundle[0]
		if left_read.aligned:
			coverage[ left_read.iv ] += 1
		if right_read.aligned:
			coverage[ right_read.iv] += 1

	"""
	Count paired-end reads
	"""

	alignment_file = HTSeq.SAM_Reader('')
	counts = collections.Counter()
	for bundle in HTSeq.pair_SAM_alignments(samfile, bundle = True):
		# Skip multiple alignments
		if len(bundle) != 1:
			continue 
		left_read, right_read = bundle[0]

		# Check to make sure the paired reads actually aligned
		if not left_read.aligned and right_read.aligned:
			count[ '_unmapped' ] += 1
			continue
		gene_ids = set()

		# Check to see if the paired reads map to the mouse genome
		for iv, val in exons[ left_read.iv ].steps():
			gene_ids |= val

		for iv, val in exons[ right_read.iv ].steps():
			gene_ids |= val

		# Check to see if the paired reads map to the spike ins 
		for iv, val in spikeins[ left_read.iv ].steps():
			gene_ids |= val

		for iv, val in spikeins[ right_read.iv ].steps():
			gene_ids |= val

		if len(gene_ids) == 1:
			gene_id = list(gene_ids)[0]
			counts[ gene_id ] += 1

		elif len(gene_ids) == 0:
			counts[ '_no_feature' ] += 1

		else:
			counts[ '_ambiguous' ] += 1 

	""" 
	Count unmapped reads
	"""
	no_unmapped_reads = counts[ '_unmapped' ]

	"""
	Count the reads mapped to the mouse genome
	"""
	no_mouse_reads = 0

	for gene_id in mouse_exon_list:
		no_mouse_reads += counts[ gene_id ]

	"""
	Count the reads mapped to the spike ins
	"""
	no_spikein_reads = 0

	for gene_id in spikein_transcript_list:
		no_spikein_reads += counts[ gene_id ]

	return coverage, counts, no_unmapped_reads, no_mouse_reads, no_spikein_reads

def load_matfiles(input_direc = None):
	file_list = os.listdir(input_direc)
	matfiles = []
	for matfile in file_list:
		if fnmatch.fnmatch(matfile, r'*.mat'):
			time = int(matfile.split('_')[0])
			dynamics_file = sio.loadmat(input_direc + matfile)
			matfiles += [[dynamics_file, time]]
	return matfiles

class dynamics_class():
	def __init__(self, matfiles):
		dynamics_dict = {}
		chip_no_dict = {}
		capture_site_dict = {} 
		num_of_cells_in_well_dict = {}
		time_point_dict = {}
		for matfile in matfiles:
			for row in xrange(matfile[0].shape[0]):
				library_id = matfile[0][row, 1]
				cell_id = matfile[0][row, 2]

				dict_loc = str(library_id) + '-' + str(cell_id)

				dynamics_dict[dict_loc] = matfile[0][row, 5:]
				chip_no_dict[dict_loc] = matfile[0][row, 0]
				capture_site_dict[dict_loc] = matfile[0][row, 3]
				num_of_cells_in_well_dict[dict_loc] = matfile[0][row, 4]
				time_point_dict[dict_loc] = matfile[1]

		self.dynamics_dict = dynamics_dict
		self.chip_no_dict = chip_no_dict
		self.capture_site_dict = capture_site_dict
		self.num_of_cells_in_well_dict = num_of_cells_in_well_dict
		self.time_point_dict = time_point_dict

class cell_object():
	def __init__(self, direc = None, samfile_name = None, dictionary = None, mouse_gtf_loc = None, spikein_gtf_loc = None):
		cell_id = parse_filename(samfile_name)
		self.NFkB_dynamics = dictionary.dynamics_dict[cell_id]
		self.chip_number = dictionary.chip_no_dict[cell_id]
		self.capture_site = dictionary.capture_site_dict[cell_id]
		self.time_point = time_point_dict[cell_id]
		self.number_of_cells = num_of_cells_in_well_dict[cell_id]

		coverage, counts, num_unmapped_reads, num_mouse_reads, num_spikein_reads = load_sequence_counts(samfile_name = None, mouse_gtf = mouse_gtf_loc, spikein_gtf = spikein_gtf_loc)
		self.coverage = coverage
		self.counts = counts
		self.num_unmapped_reads = num_unmapped_reads
		self.num_spikein_reads = num_spikein_reads
