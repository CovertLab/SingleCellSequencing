from fireworks import FireTaskBase, explicit_serialize
import seq_functions
import os

@explicit_serialize
class Count_rsem_Task(FireTaskBase):

	_fw_name = "Count_rsem_task"
	required_params = ["library_path", "aligned_name", "quant_name", "counted_name", "spikeids"]
	optional_params = []

	def run_task(self, fw_spec):

		seq_functions.count_rsem_files(self["library_path"], bammed_direc = self["aligned_name"], quant_direc = self["quant_name"], counted_direc = self["counted_name"], spikeids = self["spikeids"])

