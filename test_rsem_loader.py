import seq_functions
import os
# rsem_direc = '/scratch/PI/mcovert/dvanva/sequencing/library1/' #quant_rsem/1-10_S14.isoforms.results'
# seq_functions.load_sequence_counts_rsem(rsem_file = rsem_name, spikeids = ['AM1780SpikeIn1', 'AM1780SpikeIn4', 'AM1780SpikeIn7'])
# seq_functions.count_rsem_files(rsem_direc, spikeids = ['AM1780SpikeIn1', 'AM1780SpikeIn4', 'AM1780SpikeIn7'])


trim_direc ='/scratch/PI/mcovert/dvanva/sequencing/library1/trimmed'
file_list = os.listdir(trim_direc)

file_list_unique = []
for seq_file in file_list:
	sft = seq_file.split('.')[0]
	sft2 = sft.split('_')
	seq_file_name = sft2[0] + '_' + sft2[1]
	file_list_unique += [seq_file_name]

file_list_unique = list(set(file_list_unique))

for seq_file in file_list_unique:
	print r'Aligning ' + seq_file