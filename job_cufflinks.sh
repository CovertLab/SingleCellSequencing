# Merge stdin and stderr
#PBS -j oe

# Set maximum wall time (your job gets killed if it runs longer)
#PBS -l walltime=1:00:00

# Set maximum amount of memory your job will use (it won't get more)
#PBS -l mem=2G

# Print out useful (useless?) job-related information
echo ------------------------------------------------------
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo PBS: PYTHONPATH = $PYTHONPATH
echo ------------------------------------------------------

#########################
## Import utility functions

source "${PBS_O_WORKDIR}/util.sh"

#########################

#########################
## Check input

# Make sure the user passes a valid TOPHAT_OUTPUT_DIR environmental variable
if [ ! -d "${PBS_O_WORKDIR}/${TOPHAT_OUTPUT_DIR}" ]; then
	to_stderr "You must define the TOPHAT_OUTPUT_DIR environmental variable."
	to_stderr "And TOPHAT_OUTPUT_DIR must exist as a directory."
	exit 1
fi

# Make sure the user passes a valid CUFFLINKS_OUTPUT_DIR environmental variable
if [ ! -d "${PBS_O_WORKDIR}/${CUFFLINKS_OUTPUT_DIR}" ]; then
	to_stderr "You must define the CUFFLINKS_OUTPUT_DIR environmental variable."
	to_stderr "And CUFFLINKS_OUTPUT_DIR must exist as a directory."
	exit 1
fi

# Make sure this job knows its work allocation
if [ -z "$ARRAY_ID" ]; then
	to_stderr "ARRAY_ID environmental variable must be set"
	exit 1
fi

#########################

#########################
## Detect and set up environment on compute node

if [ -d "/state/partition1" ]; then
	WORK_DIR="/state/partition1"
else
	WORK_DIR="/tmp"
fi

TIME="$(date "+%Y%m%d.%H%M%S.%N")"

WORK_DIR="${WORK_DIR}/${TIME}.${PBS_JOBID}.${ARRAY_ID}"

CODE_DIR="$PBS_O_WORKDIR" # Assumes job submission from SingleCellSequencing

cd ${CODE_DIR}

# Get file to work on
LIST_OF_FILES="$(find "${TOPHAT_OUTPUT_DIR}" -type f -name "*.bam" -print0 | xargs -0 echo | sort)"
FULL_FILE="$(IFS=" "; set - ${LIST_OF_FILES}; shift $((ARRAY_ID - 1)); echo $1)"

# Get just the filename (strip off the directory)
BASE_FILE="$(basename $FULL_FILE)"

# Strip off the file extension
# See http://stackoverflow.com/a/965069/1647819
FILE_NO_EXT="${BASE_FILE%.*}"

# Define output directory for this cufflinks run
THIS_OUTPUT_DIR_LOCAL="${WORK_DIR}/${FILE_NO_EXT}"

#########################

stagein()
{
	echo
	echo "$(date): Copying files to work directory ${WORK_DIR}"

	mkdir -p ${WORK_DIR}

	cd ${WORK_DIR}
	scp -r "${CODE_DIR}/${FULL_FILE}" .
}

runprogram()
{
	echo "+++++ $(date): Processing ${BASE_FILE}"

	cd ${WORK_DIR}
	PATH="${PATH}:/usr/local/bin" cufflinks -o "${THIS_OUTPUT_DIR_LOCAL}" "${BASE_FILE}"
}

stageout()
{
	echo "+++++ $(date): Transferring files back"

	cd ${WORK_DIR}
	scp -r "${THIS_OUTPUT_DIR_LOCAL}" "${CODE_DIR}/${CUFFLINKS_OUTPUT_DIR}"

	echo "+++++ $(date): Cleaning up"
	echo "+++++ $(date): Removing ${BASE_FILE}"
	rm -fr "${BASE_FILE}"
	echo "+++++ $(date): Removing ${THIS_OUTPUT_DIR_LOCAL}"
	rm -fr "${THIS_OUTPUT_DIR_LOCAL}"

	cd /
	echo "+++++ $(date): Removing ${WORK_DIR}"
	rm -fr "${WORK_DIR}"
}

early()
{
	echo
	echo "##### WARNING: EARLY TERMINATION #####"
	echo
}

trap "early; stageout" 2 9 15

stagein
runprogram
stageout

exit 0