from fireworks import Firework, LaunchPad, Workflow, ScriptTask
from firetasks import Align_star_Task, Count_rsem_Task
import seq_functions
import argparse
import os
import collections
import yaml

def main(sequencing_directory, library_prefix, num_libraries, raw_data_dir):
	lpad = LaunchPad(**yaml.load(open("my_launchpad.yaml")))
	workflow_fireworks = []
	workflow_dependencies = collections.defaultdict(list)

	library_dirs = [os.path.join(sequencing_directory, library_prefix + str(i + 1)) for i in xrange(num_libraries)]
	subdirs = ["aligned_star", "quant_rsem", "counted_rsem"]

	for library_dir in library_dirs:
		seq_functions.make_directories(library_dir, subdirs)

	for library_dir in library_dirs:
		seq_functions.make_directories(library_dir, subdirs)

		name = "AlignSTAR_%s" % os.path.basename(library_dir)
		fw_align = Firework(
			[
				Align_star_Task(library_path = library_dir, trimmed_name = "trimmed", aligned_name = "aligned_star/", quant_name = "quant_rsem/")
			],
			name = name,
			spec = {"_queueadapter": {"job_name": name, "ntasks_per_node": 8, "walltime": '24:00:00'}},
			)
		workflow_fireworks.append(fw_align)

		name = "Count_%s" % os.path.basename(library_dir)
		fw_count = Firework(
			[
				Count_rsem_Task(library_path = library_dir, aligned_name = "aligned_star", quant_name = "quant_rsem", counted_name = "counted_rsem", spikeids = ['AM1780SpikeIn1', 'AM1780SpikeIn4', 'AM1780SpikeIn7'])
			],
			name = name,
			spec = {"_queueadapter": {"job_name": name}},
			)
		workflow_fireworks.append(fw_count)
		workflow_dependencies[fw_align].append(fw_count)

	lpad.add_wf(
		Workflow(workflow_fireworks, links_dict = workflow_dependencies)
		)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("sequencing_directory", help = "Directory to operate on", type = str)
	parser.add_argument("--library_prefix", help = "Prefix for library subdirectories", type = str, default = "library")
	parser.add_argument("--num_libraries", help = "Number of libraries to process", type = int, default = 10)
	parser.add_argument("--raw_data_dir", help = "Raw data directory name", type = str, default = "Raw_Data")

	args = parser.parse_args().__dict__

	main(args["sequencing_directory"], args["library_prefix"], args["num_libraries"], args["raw_data_dir"])