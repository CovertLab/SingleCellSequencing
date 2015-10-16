from fireworks import Firework, LaunchPad, Workflow, ScriptTask
from firetasks import TrimTask, AlignTask, SortTask, GeneNameTask
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
	subdirs = ['unzipped', 'trimmed', 'aligned', 'bammed', 'sorted', 'pythonized']

	for library_dir in library_dirs:
		seq_functions.make_directories(library_dir, subdirs)

		name = "Sort_%s" % os.path.basename(library_dir)
		fw_sort = Firework(
			[
				GeneNameTask(library_path = library_dir, counted_name = "counted")
			],
			name = name,
			spec = {"_queueadapter": {"job_name": name}},
			)
		workflow_fireworks.append(fw_sort)

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