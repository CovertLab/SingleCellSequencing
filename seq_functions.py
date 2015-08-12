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

def make_directories(direc_name):
	direcs_to_make = ['unzipped', 'trimmed','aligned','pythonized']
	for direc in direcs_to_make:
		if not os.path.exists(direc_name + direc):
			os.makedirs(direc_name + direc)

def unzip_file(filename, input_direc = None, output_direc = None):
	cmd = ['gunzip', input_direc + filename]
	run_cmd(cmd)
	cmd = ['mv', input_direc + filename[:-3], output_direc + filename[:-3]]
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

def load_bbmap():
	cmd = ['module', 'load', 'bbmap']
	run_cmd(cmd)

def load_STAR():
	cmd = ['module', 'load', 'STAR']
	run_cmd(cmd)

def trim_reads(direc_name, input_filename):
	left_file = direc_name + 'unzipped/' + input_filename + '_R1_001.fastq'
	right_file = direc_name + 'unzipped/' + input_filename + '_R2_001.fastq'

	ofn_temp = input_filename.split('_')
	output_filename = ofn_temp[0] + '_' + ofn_temp[1]

	output_left = direc_name + 'trimmed/' + output_filename + '_R1_trimmed.fq'
	output_right = direc_name + 'trimmed/' + output_filename + '_R2_trimmed.fq'
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

	return

def run_STAR(filename):
	option_name = []
	option = []

	option_name[0] = ' --runThreadN '
	option[0] = '8'

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
