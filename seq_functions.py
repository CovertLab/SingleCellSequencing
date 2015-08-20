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
import pyensembl
import h5py
import pandas as pd
import numpy as np


"""
Command for installing the mouse genome into pyensembl
pyensembl install --reference-name "GRCm38" --gtf-path-or-url ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz --transcript-fasta-path-or-url ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz --protein-fasta-path-or-url ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz

Load the genome using the command:

data = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')

"""

"""
Define functions
"""
def relu(x):
	return x*(x>0)

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

def sort_sam_files(direc_name, aligned_direc = 'aligned', bammed_direc = 'bammed', sorted_direc = 'sorted'):
	file_list = os.listdir(os.path.join(direc_name,aligned_direc))
	for seq_file in file_list:
		if fnmatch.fnmatch(seq_file, r'*.sam'):
			input_file = os.path.join(direc_name, aligned_direc, seq_file)

			seq_file_split = seq_file.split('_')
			reg_name = seq_file_split[0] + '_' + seq_file_split[1][:-4] + '.bam'
			sorted_name_loc = seq_file_split[0] + '_' + seq_file_split[1][:-4] + '_sorted_by_location'
			sorted_name_name = seq_file_split[0] + '_' + seq_file_split[1][:-4] + '_sorted_by_name'

			output_file = os.path.join(direc_name, bammed_direc, reg_name)
			sorted_file_loc = os.path.join(direc_name, bammed_direc, sorted_name_loc)
			sorted_file_name = os.path.join(direc_name, sorted_direc, sorted_name_name)

			cmd = ['samtools', 'view', '-bS', input_file] 
			print cmd
			bam_out = run_cmd(cmd)
			write_file(output_file, bam_out)

			#Sort by location and save
			cmd = ['samtools','sort', output_file, sorted_file_loc]
			print cmd
			run_cmd(cmd)

			#Index the bam file so it can be queried with idx stats
			cmd = ['samtools','index', sorted_file_loc + '.bam']
			print cmd
			run_cmd(cmd)

			#Sort the bam file by name and save
			cmd = ['samtools','sort', '-n', sorted_file_loc + '.bam', sorted_file_name]
			print cmd
			run_cmd(cmd)

def count_hdf5_files(direc_name, bammed_direc = 'bammed', aligned_direc = 'aligned', counted_direc = 'counted', spikeids = []):
	bammed_path = os.path.join(direc_name, bammed_direc)
	aligned_path = os.path.join(direc_name, aligned_direc)
	counted_path = os.path.join(direc_name, counted_direc)
	
	file_list = os.listdir(aligned_path)
	for seq_file in file_list:
		print seq_file, fnmatch.fnmatch(seq_file, r'*.h5')
		if fnmatch.fnmatch(seq_file, r'*.h5'):
			bfs = seq_file.split('.')
			bam_file_sorted_by_loc = bfs[0] + '_sorted_by_location.bam'
			num_mapped, num_unmapped = bam_read_count(os.path.join(bammed_path, bam_file_sorted_by_loc))
			transcripts, spikeins = load_sequence_counts_kallisto(h5file_name = os.path.join(aligned_path,seq_file), spikeids = spikeids)

			filename_save = os.path.join(counted_path, seq_file[:-3] + '.h5')
			print 'Saving ' + filename_save
			
			store = pd.HDFStore(filename_save)

			d = {'num_mapped': num_mapped, 'num_unmapped': num_unmapped}
			quality_control = pd.DataFrame(d, index = [0])

			store['quality_control'] = quality_control 
			store['transcripts'] = transcripts
			store['spikeins'] = spikeins
			store.close()

			# np.savez(filename_save, num_mapped = num_mapped, num_unmapped = num_unmapped, transcripts = transcripts, spikeins = spikeins)

def bam_read_count(bamfile_name = None):
	# Returns a tuple of the number of mapped and unmapped reads in a bam file
	
	cmd = ['samtools', 'idxstats', bamfile_name]
	print cmd
	# To read the output of idxstats line by line, we need to use Popen without communicate()
	environ = {
		"PATH": os.environ["PATH"],
		"LANG": "C",
		"LC_ALL": "C",
		}
	output = subprocess.Popen(cmd, stdout = subprocess.PIPE, env = environ)

	mapped = 0
	unmapped = 0

	for line in output.stdout:
		rname, rlen, nm, nu = line.rstrip().split()
		mapped += relu(int(nm) - int(nu))
		if nm > 0:
			unmapped += 2* int(nu)
		else:
			unmapped += int(nu)
	print mapped, unmapped
	return (mapped, unmapped)

def load_ensembl_gene_ids(mouse_gtf = None):
	data = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')
	return data

def load_sequence_counts_kallisto(h5file_name = None, spikeids = []):
	"""
	To make: dictionary with counts, addressable by gene id name
	dictionary with bootstrap errors, addressable by gene id name
	First load HDF file
	Create dictionary 
	"""

	# Load mouse genome with pyensembl
	# mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')

	seq_file = h5py.File(h5file_name, 'r')
	class Kdata(object):
		def __init__(self,f):
			self.kversion = ''
			self.idxversion = 0
			self.num_bootstraps = 0
			self.start_time = ''
			self.call = ''

			num_targets = len(f['aux']['ids'])
			self.ids = np.array(f['aux']['ids'])

			# self.gene_names = np.zeros(len(self.ids), dtype = str)
			# for i in xrange(len(self.ids)):
			# 	print i
			# 	if self.ids[i] not in spikeids:
			# 		self.gene_names[i] = mouse_genome.gene_name_of_transcript_id(self.ids[i])
			
			self.lengths = np.array(f['aux']['lengths'])
			self.eff_lengths = np.array(f['aux']['eff_lengths'])

			return None

	def get_aux_data(f):
		kdata = Kdata(f)
		kdata.kversion = f['aux']['kallisto_version'][0]
		kdata.idxversion = f['aux']['index_version'][0]
		kdata.num_targets = len(f['aux']['ids'])
		kdata.num_bootstrap = f['aux']['num_bootstrap'][0]
		kdata.start_time = f['aux']['start_time'][0]
		kdata.call = f['aux']['call'][0]
		return kdata

	bs_col_list = []
	kdata = get_aux_data(seq_file)

	d = {'target_id': kdata.ids, 'length': kdata.lengths, 'eff_length': kdata.eff_lengths}
	data_frame = pd.DataFrame(d)
	est_counts = pd.Series(seq_file['est_counts'], index = data_frame.index)
	data_frame['est_counts'] = est_counts

	if kdata.num_bootstrap > 0:
		# Add the boot strap count estimates as columns
		for i in xrange(kdata.num_bootstrap):
			bs_key = 'b' + str(i)
			bs_col_list.append(bs_key)
			bs_counts = pd.Series(f['bootstrap']['bs_key'], index = data_frame.index)
			data_frame[bs_key] = bs_counts

		# Compute statistics
		data_frame['est_counts_mean'] = data_frame[bs_col_list].mean(axis = 1)
		data_frame['est_counts_stdev'] = data_frame[bs_col_list].std(axis = 1)
		data_frame['est_counts_sem'] = data_frame[bs_col_list].sem(axis = 1)
		data_frame['est_counts_min'] = data_frame[bs_col_list].min(axis = 1)
		data_frame['est_counts_med'] = data_frame[bs_col_list].median(axis = 1)
		data_frame['est_counts_max'] = data_frame[bs_col_list].max(axis = 1)

	spike_ins = data_frame[data_frame['target_id'].isin(spikeids)]
	spike_locs = []
	for spikeid in spikeids:
		spike_locs += np.where(kdata.ids == spikeid)
	transcript_ids_only = np.delete(kdata.ids,spike_locs)
	transcripts = data_frame[data_frame['target_id'].isin(transcript_ids_only)]

	# Compute TPM
	if kdata.num_bootstrap == 0:
		divisor = (transcripts.est_counts / transcripts.eff_length).sum()
		transcripts['tpm'] = ( (transcripts.est_counts / transcripts.eff_length) / divisor) * 1e6
	else:
		divisor = (transcripts.est_counts_mean / transcripts.eff_length).sum()
		transcripts['tpm'] = ( (transcripts.est_counts_mean / transcripts.eff_length) / divisor) * 1e6

	transcripts.set_index('target_id', inplace = True)
	spike_ins.set_index('target_id', inplace = True)

	return transcripts, spike_ins

def load_sequence_counts_STAR(bamfile_name = None, mouse_gtf = '/scratch/PI/mcovert/dvanva/sequencing/ref_seq_cdna/Mus_musculus.GRCm38.81.gtf', spikein_gtf = '/scratch/PI/mcovert/dvanva/sequencing/ref_seq/spikeInsAM1780.gtf'):
		
	"""
	Load gtf files
	"""

	mouse_gtf_file = HTSeq.GFF_Reader(mouse_gtf)
	spikein_gtf_file = HTSeq.GFF_Reader(spikein_gtf)

	"""
	Load SAM file
	"""
	samfile = HTSeq.BAM_Reader(bamfile_name)

	"""
	Construct list of genes for the mouse genome and spikeins
	"""
	
	spikein_transcript_list = []
	for feature in spikein_gtf_file:
		if feature.type == 'exon':
			spikein_transcript_list += [feature.attr['gene_id']]

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
	num_unmapped_reads = counts[ '_unmapped' ]

	"""
	Count the reads mapped to the exons in the mouse genome
	"""
	num_mouse_reads = 0

	for gene_id in mouse_exon_list:
		num_mouse_reads += counts[ gene_id ]

	"""
	Count the reads mapped to the spike ins
	"""
	num_spikein_reads = 0

	for gene_id in spikein_transcript_list:
		num_spikein_reads += counts[ gene_id ]

	return coverage, counts, num_mouse_reads, num_spikein_reads

def load_matfiles(input_direc = None):
	file_list = os.listdir(input_direc)
	matfiles = []
	var_list = ['All300minsdata', 'All150minsNoStimulationdata', 'All150minsdata','All75minsdata','All0minsdata'] 
	for it in xrange(len(file_list)):
		matfile = file_list[it]
		if fnmatch.fnmatch(matfile, r'*.mat'):
			time = int(matfile.split('_')[0])
			dynamics_file = sio.loadmat(os.path.join(input_direc,matfile))
			if var_list[it] == 'All150minsNoStimulationdata':
				condition = 'NoStim'
			else:
				condition = 'Stim'
			matfiles += [[dynamics_file[var_list[it]], time, condition]]
	return matfiles

class dynamics_class():
	def __init__(self, matfiles):
		dynamics_dict = {}
		chip_no_dict = {}
		capture_site_dict = {} 
		num_of_cells_in_well_dict = {}
		time_point_dict = {}
		condition_dict = {}
		for matfile in matfiles:
			for row in xrange(matfile[0].shape[0]):
				library_id = int(matfile[0][row, 1])
				cell_id = int(matfile[0][row, 2])

				dict_loc = str(library_id) + '-' + str(cell_id)

				dynamics_dict[dict_loc] = matfile[0][row, 5:]
				condition_dict[dict_loc] = matfile[2]
				chip_no_dict[dict_loc] = matfile[0][row, 0]
				capture_site_dict[dict_loc] = matfile[0][row, 3]
				num_of_cells_in_well_dict[dict_loc] = matfile[0][row, 4]
				time_point_dict[dict_loc] = matfile[1]

		self.dynamics_dict = dynamics_dict
		self.chip_no_dict = chip_no_dict
		self.capture_site_dict = capture_site_dict
		self.num_of_cells_in_well_dict = num_of_cells_in_well_dict
		self.time_point_dict = time_point_dict
		self.condition_dict = condition_dict

class cell_object():
	def __init__(self, h5_file = None, dictionary = None):
		cell_id = parse_filename(os.path.basename(h5_file))
		self.NFkB_dynamics = dictionary.dynamics_dict[cell_id]
		self.chip_number = dictionary.chip_no_dict[cell_id]
		self.capture_site = dictionary.capture_site_dict[cell_id]
		self.time_point = dictionary.time_point_dict[cell_id]
		self.number_of_cells = dictionary.num_of_cells_in_well_dict[cell_id]
		self.condition = dictionary.condition_dict[cell_id]
		self.clusterID = 0
		self.quality = 0

		seq_data = pd.HDFStore(h5_file)
		self.num_mapped = seq_data['quality_control']['num_mapped']
		self.num_unmapped = seq_data['quality_control']['num_unmapped']
		self.transcripts = seq_data['transcripts']
		self.spikeins = seq_data['spikeins']
		self.total_estimated_counts = self.transcripts.est_counts.sum() + self.spikeins.est_counts.sum()
		seq_data.close()


