from fireworks import FireTaskBase, explicit_serialize
import seq_functions
import os

@explicit_serialize
class SortTask(FireTaskBase):

	_fw_name = "SortTask"
	required_params = ["library_path", "aligned_name", "bammed_name", "sorted_name"]
	optional_params = []

	def run_task(self, fw_spec):

		seq_functions.sort_sam_files(self["library_path"], aligned_direc = self["aligned_name"], bammed_direc = self["bammed_name"], sorted_direc = self["sorted_name"])