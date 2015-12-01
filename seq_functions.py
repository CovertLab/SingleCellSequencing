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
import scipy.cluster.hierarchy as sch
import rpy2
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
rpy2.robjects.numpy2ri.activate()


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

def count_reads(file_name):
	cmd = ['wc', 'l', file_name]
	counts = int(run_cmd(cmd).split()[0])

def run_star_rsem(direc_name, input_filename, genomeDir = '/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna', rsem_reference_name = '/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/rsem/rsem_mouse_genome_ref', trimmed_name = 'trimmed', aligned_name = 'aligned_star/', quant_name = 'quant_rsem/'):
	
	# Align with STAR and quantify with rsem
	input_directory = os.path.join(direc_name,trimmed_name)
	aligned_directory = os.path.join(direc_name, aligned_name)
	quant_directory = os.path.join(direc_name, quant_name)

	pair1 = os.path.join(input_directory, input_filename+'_R1_trimmed.fq')
	pair2 = os.path.join(input_directory, input_filename+'_R2_trimmed.fq')

	cmd = ['rsem-calculate-expression', '--star', '--keep-intermediate-files', '--paired-end', pair1, pair2, rsem_reference_name, os.path.join(aligned_directory, input_filename)]
	print cmd
	run_cmd(cmd)

	genes_name = os.path.join(aligned_directory, input_filename +'.genes.results')
	genes_moved = os.path.join(quant_directory, input_filename +'.genes.results')

	cmd = ['mv', genes_name, genes_moved]
	print cmd
	run_cmd(cmd)

	isoform_name = os.path.join(aligned_directory, input_filename +'.isoforms.results')
	isoform_moved = os.path.join(quant_directory, input_filename +'.isoforms.results')

	cmd = ['mv', isoform_name, isoform_moved]
	print cmd
	run_cmd(cmd)


def generate_genome_index():
	option_name = []
	option = []

	option_name += option_name + ['--runThreadN']
	option += ['16']
	option_name += ['--runMode'] 
	option += ['genomeGenerate']
	option_name += ['--genomeDir']
	option += ['/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna']
	option_name += ['--genomeFastaFiles']
	option += ['/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/mouse_genome_w_spikeins.fa']
	option_name += ['--sjdbGTFfile']
	option += ['/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/mouse_genome_w_spikeins.gtf']
	option_name += ['--sjdbOverhang ']
	option += ['76']

	cmd = ['STAR']
	for j in xrange(len(option_name)):
		cmd = cmd + [option_name[j]] + [option[j]]

	print cmd
	run_cmd(cmd)


def run_STAR_direc(direc_name, trimmed_name = 'trimmed', aligned_name = 'aligned_star'):
	output_directory = os.path.join(direc_name, aligned_name)
	file_list = os.listdir(os.path.join(direc_name, trimmed_name))
	for seq_file in file_list:
		sft = seq_file.split('.')[0]
		sft2 = sft.split('_')
		seq_file_name = sft[0] + '_' + sft2[1]
		if not os.path.exists(os.path.join(output_directory, seq_file_name + '.sam')):
			run_STAR(direc_name, seq_file_name)

def generate_rsem_reference_genome():
	cmd = ['rsem-prepare-reference','--star','--gtf','/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/mouse_genome_w_spikeins.gtf','/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/mouse_genome_w_spikeins.fa','/scratch/PI/mcovert/dvanva/sequencing/ref_seq_dna/rsem/rsem_mouse_genome_ref']
	print cmd
	run_cmd(cmd)


def run_kallisto_index(input_filename, output_filename):
	cmd = ['kallisto', 'index', '-i', output_filename, input_filename]
	run_cmd(cmd)
	for seq_file in file_list:
		sft = seq_file.split('.')[0]
		sft2 = sft.split('_')
		seq_file_name = sft[0] + '_' + sft2[1]
		if not os.path.exists(os.path.join(output_directory, seq_file_name + '.sam')):
			run_STAR(direc_name, seq_file_name)

def kallisto_direc(direc_name):
	output_directory = os.path.join(direc_name, 'aligned')
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

	mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')
	
	for seq_file_temp in file_list:
		if fnmatch.fnmatch(seq_file_temp, r'*.h5'):
			seq_file = seq_file_temp
	transcripts, spikeins = load_sequence_counts_kallisto(h5file_name = os.path.join(aligned_path,seq_file), spikeids = spikeids)

	transcript_names = transcripts.index
	gene_names = []
	counter = 0
	for transcript in transcript_names:
		counter += 1
		print counter
		gene_names += [mouse_genome.gene_name_of_transcript_id(transcript)]
	unique_gene_names = list(set(gene_names))

	transcript_dict = {}
	for gene in unique_gene_names:
		transcript_dict[gene] = []

	for transcript in transcript_names:
		gene_name = mouse_genome.gene_name_of_transcript_id(transcript)
		transcript_dict[gene_name] += [transcript]

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

			# Create a transcript table indexed by gene name
			counts_list = []
			fpkm_list = []
			tpm_list = []
			length_list = []
			for j in xrange(len(unique_gene_names)):
				gene = unique_gene_names[j]
				transcript_list = transcript_dict[gene]
				counts = 0
				fpkm = 0
				tpm = 0
				mean_eff_length = 0
				for transcript in transcript_list:
					counts += transcripts.loc[transcript]['est_counts']
					fpkm += transcripts.loc[transcript]['fpkm']
					tpm += transcripts.loc[transcript]['tpm']
					mean_eff_length += transcripts.loc[transcript]['eff_length'] / len(transcript_list)
				counts_list += [counts]
				fpkm_list += [fpkm]
				tpm_list += [tpm]
				length_list += [mean_eff_length]

			d = {'gene_name': unique_gene_names, 'est_counts': counts_list, 'fpkm': fpkm_list, 'tpm': tpm_list, 'mean_eff_length': length_list}
			gene_counts = pd.DataFrame(d)
			gene_counts.set_index('gene_name', inplace = True)

			# Store data frames in HDF5 format
			store['quality_control'] = quality_control 
			store['transcripts'] = transcripts
			store['gene_counts'] = gene_counts
			store['spikeins'] = spikeins
			store.close()

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
	transcripts = data_frame.loc[data_frame['target_id'].isin(transcript_ids_only)]
	transcripts.is_copy = False

	# Compute TPM
	if kdata.num_bootstrap == 0:
		divisor = (transcripts['est_counts'] / transcripts['eff_length']).sum()
		tpm = ((transcripts['est_counts'] / transcripts['eff_length']) / divisor) * 1e6
		fpkm = transcripts['est_counts'] / (transcripts['eff_length'] * transcripts['est_counts'].sum()) * 1e9
		transcripts['fpkm'] = fpkm
		transcripts['tpm'] = tpm

	else:
		divisor = (transcripts['est_counts_mean'] / transcripts['eff_length']).sum()
		tpm = ( (transcripts['est_counts_mean'] / transcripts['eff_length']) / divisor) * 1e6
		fpkm = transcripts['est_counts_mean'] / (transcripts['eff_length'] * transcripts['est_counts_mean'].sum()) * 1e9
		transcripts['fpkm'] = fpkm
		transcripts['tpm'] = tpm

	transcripts.set_index('target_id', inplace = True)
	spike_ins.set_index('target_id', inplace = True)

	return transcripts, spike_ins

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

def quality_control(list_of_cells, min_num_of_reads = 200000, max_fraction_spikein = 0.3, max_fraction_unmapped = 0.5):
	for cell in list_of_cells:
		cell.quality = 1

		# Screen the number of mapped reads
		if cell.num_mapped < min_num_of_reads:
			cell.quality = 0

		# Screen the fraction of reads that map to spike ins
		spikein_counts = cell.spikeins.loc['Spike1']['est_counts'] + cell.spikeins.loc['Spike4']['est_counts'] + cell.spikeins.loc['Spike7']['est_counts']
		if spikein_counts/cell.total_estimated_counts > max_fraction_spikein:
			cell.quality = 0

		# Remove cells where 50% or more of the reads were unmapped
		if cell.num_unmapped/cell.num_mapped > max_fraction_unmapped:
			cell.quality = 0

	new_list_of_cells = []
	for cell in list_of_cells:
		if cell.quality == 1:
			new_list_of_cells += [cell]

	return new_list_of_cells

def cell_cluster(list_of_cells, max_clusters = 3):
	R = rpy2.robjects.r
	DTW = importr('dtw')
	longest_time = 0
	number_of_cells = 0
	for cell in list_of_cells:
		number_of_cells += 1
		longest_time = np.amax([longest_time, cell.NFkB_dynamics.shape[0]])

	dynamics_matrix = np.zeros((number_of_cells,longest_time))

	"""
	Fill up the heat map matrix
	"""

	cell_counter = 0
	for j in xrange(number_of_cells):
		cell = list_of_cells[j]		
		dynam = cell.NFkB_dynamics
		dynamics_matrix[cell_counter,0:dynam.shape[0]] = dynam
		cell_counter += 1

	"""
	Perform hierarchical clustering
	"""

	distance_matrix = np.zeros((number_of_cells, number_of_cells))
	for i in xrange(number_of_cells):
		for j in xrange(number_of_cells):
			alignment = R.dtw(dynamics_matrix[i,:], dynamics_matrix[j,:], keep = True)
			distance_matrix[i,j] = alignment.rx('distance')[0][0]

	Y = sch.linkage(distance_matrix, method = 'centroid')

	clusters = sch.fcluster(Y, max_clusters, criterion = 'maxclust')

	for j in xrange(number_of_cells):
		list_of_cells[j].clusterID = clusters[j]

	return list_of_cells

def remove_unidentified_genes(list_of_cells):
	cell = list_of_cells[0]
	gene_keys = cell.fpkm.index
	print "Removing unidentified genes ..."
	
	counter = 0

	zero_genes = set(cell.fpkm[cell.fpkm == 0].index.tolist())
	for cell in list_of_cells:
		zero_genes_new = set(cell.fpkm[cell.fpkm == 0].index.tolist())
		zero_genes &= zero_genes_new 

	print str(len(zero_genes)) + ' of ' + str(len(gene_keys)) + ' genes not detected'
	for cell in list_of_cells:
		cell.fpkm.drop(list(zero_genes), axis = 0, inplace = True)

	return list_of_cells

def remove_unidentified_genes_rsem(list_of_cells):
	cell = list_of_cells[0]
	gene_keys = cell.fpkm_rsem.index
	print "Removing unidentified genes ..."
	
	counter = 0

	zero_genes = set(cell.fpkm[cell.fpkm_rsem == 0].index.tolist())
	for cell in list_of_cells:
		zero_genes_new = set(cell.fpkm[cell.fpkm_rsem == 0].index.tolist())
		zero_genes &= zero_genes_new 

	print str(len(zero_genes)) + ' of ' + str(len(gene_keys)) + ' genes not detected'
	for cell in list_of_cells:
		cell.fpkm_rsem.drop(list(zero_genes), axis = 0, inplace = True)

	return list_of_cells

def remove_jackpotting_genes(list_of_cells, jackpot_threshold = 0.05):
	print "Removing jackpotting from transcripts ..."
	cell = list_of_cells[0]
	gene_keys = cell.fpkm.index
	
	num_removed_genes = 0
	for cell in list_of_cells:

		jackpotted_genes = []
		fraction_of_transcriptome = cell.fpkm/cell.fpkm.sum()
		jackpotted_genes = fraction_of_transcriptome[fraction_of_transcriptome > jackpot_threshold].index.tolist()

		num_removed_genes += len(jackpotted_genes)
		total_fpkm = cell.fpkm.sum()
		for gene in jackpotted_genes:
			total_fpkm -= cell.fpkm.loc[gene]
		mean_fpkm = total_fpkm/(len(gene_keys) - len(jackpotted_genes))
		for gene in jackpotted_genes:
			cell.fpkm.loc[gene] = mean_fpkm

	print str(num_removed_genes) + ' jackpotting events identified and removed'

	return list_of_cells

def add_tpm_normalization(list_of_cells):
	print "Adding tpm normalization to cell objects ..."
	for cell in list_of_cells:
		cell.tpm = cell.fpkm/cell.fpkm.sum() * 1e6
	return list_of_cells

def add_tpm_normalization_rsem(list_of_cells):
	print "Adding tpm normalization to cell objects ..."
	for cell in list_of_cells:
		cell.tpm_rsem = cell.fpkm_rsem/cell.fpkm_rsem.sum() * 1e6
	return list_of_cells

def remove_low_expression_genes(list_of_cells, tpm_plus_one_threshold = 100):
	cell = list_of_cells[0]
	gene_keys = cell.tpm.index
	num_of_cells = len(list_of_cells)

	dropout_frequency = pd.DataFrame(np.zeros((len(gene_keys),2), dtype = 'float32'), index = gene_keys, columns = ['tpm','dropout'])

	for cell in list_of_cells:
		dropout_frequency.loc[:,'tpm'] += (1+cell.tpm)/num_of_cells
		zero_genes = cell.tpm[cell.tpm == 0].index.tolist()
		dropout_frequency.ix[zero_genes,'dropout'] += np.float32(num_of_cells) ** -1

	low_expression_genes = dropout_frequency[dropout_frequency['tpm'] + 1 < 100].index.tolist()

	for cell in list_of_cells:
		cell.fpkm.drop(list(low_expression_genes), axis = 0, inplace = True)
		cell.tpm.drop(list(low_expression_genes), axis = 0, inplace = True)

	return list_of_cells
	
class dynamics_class():
	def __init__(self, matfiles):
		dynamics_dict = {}
		chip_no_dict = {}
		capture_site_dict = {} 
		num_of_cells_in_well_dict = {}
		time_point_dict = {}
		condition_dict = {}
		library_number_dict = {}
		cell_number_dict = {}
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
				cell_number_dict[dict_loc] = cell_id
				library_number_dict[dict_loc] = library_id

		self.library_number_dict = library_number_dict
		self.cell_number_dict = cell_number_dict
		self.dynamics_dict = dynamics_dict
		self.chip_no_dict = chip_no_dict
		self.capture_site_dict = capture_site_dict
		self.num_of_cells_in_well_dict = num_of_cells_in_well_dict
		self.time_point_dict = time_point_dict
		self.condition_dict = condition_dict

def count_rsem_files(direc_name, bammed_direc = 'aligned_star', quant_direc = 'quant_rsem', counted_direc = 'counted_rsem', spikeids = []):
	quant_path = os.path.join(direc_name, quant_direc)
	counted_path = os.path.join(direc_name, counted_direc)
	bammed_path = os.path.join(direc_name, bammed_direc)

	file_list = os.listdir(quant_path)

	mouse_genome = pyensembl.Genome(reference_name = 'GRCm38', gtf_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/gtf/mus_musculus/Mus_musculus.GRCm38.81.gtf.gz', transcript_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz', protein_fasta_path_or_url = 'ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/pep/Mus_musculus.GRCm38.pep.all.fa.gz')
	
	for seq_file_temp in file_list:
		if fnmatch.fnmatch(seq_file_temp, r'*.isoforms.results'):
			seq_file = seq_file_temp

	gene_names = mouse_genome.gene_names()

	transcript_dict = {}
	counter = 0
	for gene in gene_names:
		transcript_dict[gene] = mouse_genome.transcript_ids_of_gene_name(gene)
		counter += 1
		print counter

	for seq_file in file_list:
		print seq_file, fnmatch.fnmatch(seq_file, r'*.isoforms.results')
		if fnmatch.fnmatch(seq_file, r'*.isoforms.results'):
			bfs = seq_file.split('.')
			bam_file_sorted_by_loc = bfs[0] + '.transcript.sorted.bam'
			num_mapped, num_unmapped = bam_read_count(os.path.join(bammed_path, bam_file_sorted_by_loc))
			transcripts, spikeins = load_sequence_counts_rsem(rsem_file = os.path.join(quant_path,seq_file), spikeids = spikeids)

			filename_save = os.path.join(counted_path, seq_file + '.h5')
			print 'Saving ' + filename_save
			
			store = pd.HDFStore(filename_save)

			d = {'num_mapped': num_mapped, 'num_unmapped': num_unmapped}
			quality_control = pd.DataFrame(d, index = [0])

			# Create a transcript table indexed by gene name - sum over all isoforms
			counts_list = []
			fpkm_list = []
			tpm_list = []
			for j in xrange(len(gene_names)):
				gene = gene_names[j]
				transcript_list = transcript_dict[gene]
				counts = 0
				fpkm = 0
				tpm = 0
				for transcript in transcript_list:
					if transcript in transcripts.index:
						counts += transcripts.loc[transcript]['est_counts']
						fpkm += transcripts.loc[transcript]['fpkm']
						tpm += transcripts.loc[transcript]['tpm']
				counts_list += [counts]
				fpkm_list += [fpkm]
				tpm_list += [tpm]

			d = {'gene_name': gene_names, 'est_counts': counts_list, 'fpkm': fpkm_list, 'tpm': tpm_list}
			gene_counts = pd.DataFrame(d)
			gene_counts.set_index('gene_name', inplace = True)

			# Store data frames in HDF5 format
			store['quality_control'] = quality_control 
			store['transcripts'] = transcripts
			store['gene_counts'] = gene_counts
			store['spikeins'] = spikeins
			store.close()


def load_sequence_counts_rsem(rsem_file = None, spikeids = []):
	data_frame = pd.read_csv(rsem_file, sep = r'\t')
	spike_mask = data_frame['transcript_id'].isin(spikeids)
	spikeins = data_frame[spike_mask]
	spikeins.is_copy = False

	transcript_mask = np.logical_not(spike_mask)
	transcripts = data_frame[transcript_mask]
	transcripts.is_copy = False

	transcripts = transcripts[transcripts['effective_length'] > 0]

	fpkm = transcripts['expected_count'] / (transcripts['effective_length'] * transcripts['expected_count'].sum()) * 1e9
	divisor = (transcripts['expected_count'] / transcripts['effective_length']).max()
	tpm = ( (transcripts['expected_count'] / transcripts['effective_length']) / divisor) * 1e6


	transcripts['FPKM'] = fpkm
	transcripts['TPM'] = tpm

	transcripts.rename(columns = {'transcript_id' : 'target_id', 'FPKM': 'fpkm', 'TPM': 'tpm', 'effective_length': 'eff_length', 'expected_count': 'est_counts'}, inplace = True)
	spikeins.rename(columns = {'gene_id' : 'target_id', 'FPKM': 'fpkm', 'TPM': 'tpm', 'effective_length': 'eff_length', 'expected_count': 'est_counts'}, inplace = True)

	transcripts.set_index('target_id', inplace = True)
	spikeins.set_index('target_id', inplace = True)

	print transcripts
	print spikeins

	return transcripts, spikeins


class cell_object():
	def __init__(self, h5_file = None, dictionary = None):

		cell_id = parse_filename(os.path.basename(h5_file))

		# Load experimental metadata
		self.id = cell_id
		self.library_number = dictionary.library_number_dict[cell_id]
		self.cell_number = dictionary.cell_number_dict[cell_id]
		self.chip_number = dictionary.chip_no_dict[cell_id]
		self.capture_site = dictionary.capture_site_dict[cell_id]
		self.time_point = dictionary.time_point_dict[cell_id]
		self.number_of_cells = dictionary.num_of_cells_in_well_dict[cell_id]
		self.condition = dictionary.condition_dict[cell_id]
		self.clusterID = 0
		self.quality = 0

		# Load signaling dynamics
		self.NFkB_dynamics = dictionary.dynamics_dict[cell_id]

		# Load transcriptome
		seq_data = pd.HDFStore(h5_file)
		self.num_mapped = seq_data['quality_control']['num_mapped'].item()/2
		self.num_unmapped = seq_data['quality_control']['num_unmapped'].item()/2
		
		# Un comment if we want to keep isoform level information
		# self.transcripts = seq_data['transcripts']

		# Include data summed over isoforms - be wary of using the counts directly
		self.fpkm = seq_data['gene_counts'].loc[:,'fpkm']
		self.est_counts = seq_data['gene_counts'].loc[:,'est_counts']
		self.mean_eff_length = seq_data['gene_counts'].loc[:,'mean_eff_length']
		self.spikeins = seq_data['spikeins']
		self.total_estimated_counts = seq_data['transcripts'].est_counts.sum() + seq_data['spikeins'].est_counts.sum()
		seq_data.close()

class cell_object_rsem():
	def __init__(self, h5_kallisto_file = None, h5_rsem_file = None, dictionary = None):

		cell_id = parse_filename(os.path.basename(h5_kallisto_file))

		# Load experimental metadata
		self.id = cell_id
		self.library_number = dictionary.library_number_dict[cell_id]
		self.cell_number = dictionary.cell_number_dict[cell_id]
		self.chip_number = dictionary.chip_no_dict[cell_id]
		self.capture_site = dictionary.capture_site_dict[cell_id]
		self.time_point = dictionary.time_point_dict[cell_id]
		self.number_of_cells = dictionary.num_of_cells_in_well_dict[cell_id]
		self.condition = dictionary.condition_dict[cell_id]
		self.clusterID = 0
		self.quality = 0

		# Load signaling dynamics
		self.NFkB_dynamics = dictionary.dynamics_dict[cell_id]

		# Load transcriptome
		seq_data_kal = pd.HDFStore(h5_kallisto_file)
		seq_data_rsem = pd.HDFStore(h5_rsem_file)

		self.num_mapped = seq_data_kal['quality_control']['num_mapped'].item()/2
		self.num_unmapped = seq_data_kal['quality_control']['num_unmapped'].item()/2

		self.num_mapped_rsem = seq_data_rsem['quality_control']['num_mapped'].item()/2
		self.num_unmapped_rsem = seq_data_rsem['quality_control']['num_unmapped'].item()/2
		
		# Un comment if we want to keep isoform level information
		# self.transcripts = seq_data['transcripts']

		# Include data summed over isoforms - be wary of using the counts directly
		self.fpkm = seq_data_kal['gene_counts'].loc[:,'fpkm']
		self.spikeins = seq_data_kal['spikeins']
		self.total_estimated_counts = seq_data_kal['transcripts'].est_counts.sum() + seq_data_kal['spikeins'].est_counts.sum()
		
		self.fpkm_rsem = seq_data_rsem['gene_counts'].loc[:,'fpkm']
		self.spikeins_rsem = seq_data_rsem['spikeins']
		self.total_estimated_counts_rsem = seq_data_rsem['transcripts'].est_counts.sum() + seq_data_rsem['spikeins'].est_counts.sum()

		seq_data_kal.close()
		seq_data_rsem.close()


