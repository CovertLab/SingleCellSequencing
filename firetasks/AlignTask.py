from fireworks import FireTaskBase, explicit_serialize
import seq_functions
import os

@explicit_serialize
class AlignTask(FireTaskBase):

	_fw_name = "AlignTask"
	required_params = ["library_path", "trimmed_name", "aligned_name"]
	optional_params = []

	def run_task(self, fw_spec):

		trim_direc = os.path.join(self["library_path"], self["trimmed_name"])
		file_list = os.listdir(trim_direc)
		for seq_file in file_list:
			sft = seq_file.split('.')[0]
			sft2 = sft.split('_')
			seq_file_name = sft2[0] + '_' + sft2[1]
			print r'Aligning ' + seq_file_name
			seq_functions.run_kallisto_quant(self["library_path"], seq_file_name, trimmed_name = self["trimmed_name"], aligned_name = self["aligned_name"])
