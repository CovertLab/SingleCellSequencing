from fireworks import FireTaskBase, explicit_serialize
import seq_functions
import os

@explicit_serialize
class TrimTask(FireTaskBase):

	_fw_name = "TrimTask"
	required_params = ["library_path", "unzipped_name", "trimmed_name"]
	optional_params = []

	def run_task(self, fw_spec):

		unzip_direc = os.path.join(self["library_path"], self["unzipped_name"])

		file_list = os.listdir(unzip_direc)
		for seq_file in file_list:
			sft = seq_file.split('.')[0]
			sft2 = sft.split('_')
			seq_file_name = sft2[0] + '_' + sft2[1] + '_' + sft2[2] 
			print r'Trimming ' + seq_file_name
			seq_functions.trim_reads(self["library_path"], seq_file_name, unzipped_name = self["unzipped_name"], trimmed_name = self["trimmed_name"])