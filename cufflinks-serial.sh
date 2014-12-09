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

echo "$(date): Running $(basename $0)"

#########################
## Process files

echo "$(date): Processing BAM files in ${TOPHAT_OUTPUT_DIR}"
echo

echo "========================="

find "${TOPHAT_OUTPUT_DIR}" -type f -name "*.bam" -print0 | sort -z | \
while IFS='' read -r -d '' FULL_FILE; do

	# Get just the filename (strip off the directory)
	BASE_FILE="$(basename $FULL_FILE)"

	# Strip off the file extension
	# See http://stackoverflow.com/a/965069/1647819
	FILE_NO_EXT="${BASE_FILE%.*}"

	# Define output directory for this cufflinks run
	THIS_OUTPUT_DIR="${CUFFLINKS_OUTPUT_DIR}/${FILE_NO_EXT}"

	echo "+++++ $(date): Processing ${FULL_FILE}"

	cufflinks -o "${THIS_OUTPUT_DIR}" "${FULL_FILE}"

	echo
done

echo "========================="
echo

#########################

echo "$(date): Done running $(basename $0)"