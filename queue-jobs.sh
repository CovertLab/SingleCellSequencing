#!/bin/bash

#########################
## Import utility functions

source util.sh

#########################

#########################
## Check input

# Make sure the user passes a valid TOPHAT_OUTPUT_DIR environmental variable
if [ ! -d "${TOPHAT_OUTPUT_DIR}" ]; then
	to_stderr "You must define the TOPHAT_OUTPUT_DIR environmental variable."
	to_stderr "And TOPHAT_OUTPUT_DIR must exist as a directory."
	exit 1
fi

# Make sure the user passes a valid CUFFLINKS_OUTPUT_DIR environmental variable
if [ ! -d "${CUFFLINKS_OUTPUT_DIR}" ]; then
	to_stderr "You must define the CUFFLINKS_OUTPUT_DIR environmental variable."
	to_stderr "And CUFFLINKS_OUTPUT_DIR must exist as a directory."
	exit 1
fi

#########################

#########################
## Submit jobs

# Determine number of BAM files
NUMBER_BAM_FILES="$(find "${TOPHAT_OUTPUT_DIR}" -type f -name "*.bam" -printf x | wc -c)"

for (( i=1; i<=$NUMBER_BAM_FILES; i++ )); do
	THIS_CUFFLINKS_JOB=$(qsub -v TOPHAT_OUTPUT_DIR="${TOPHAT_OUTPUT_DIR}",\
CUFFLINKS_OUTPUT_DIR="${CUFFLINKS_OUTPUT_DIR}",\
ARRAY_ID="${i}" ./job_cufflinks.sh)
	echo "Submitted ${THIS_CUFFLINKS_JOB}"
done

#########################