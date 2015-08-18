from fireworks import FireTaskBase, explicit_serialize
import seq_functions
import os

@explicit_serialize
class CountTask(FireTaskBase):

	_fw_name = "CountTask"
	required_params = ["library_path", "bammed_name", "sorted_name", "counted_name"]
	optional_params = []

	def run_task(self, fw_spec):

		seq_functions.count_bam_files(self["library_path"], bammed_direc = self["bammed_name"], sorted_direc = self["sorted_name"], counted_direc = self["counted_name"])