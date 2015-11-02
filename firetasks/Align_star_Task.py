from fireworks import FireTaskBase, explicit_serialize
import seq_functions
import os

@explicit_serialize
class Align_star_Task(FireTaskBase):

	_fw_name = "Align_star_task"
	required_params = ["library_path", "trimmed_name", "aligned_name", "quant_name"]
	optional_params = []

	def run_task(self, fw_spec):

		trim_direc = os.path.join(self["library_path"], self["trimmed_name"])
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
			seq_functions.run_star_rsem(self["library_path"], seq_file, trimmed_name = self["trimmed_name"], aligned_name = self["aligned_name"], quant_name = self["quant_name"])

