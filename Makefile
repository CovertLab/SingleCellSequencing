.PHONY: cufflinks-serial queue-jobs

# Protect against environmental variable not being set
guard-env-%:
	@ if [ "${${*}}" == "" ]; then \
		echo "Environment variable $* not set"; \
		exit 1; \
	fi

# Protect against directory specified by environmental variable not existing
guard-dir-%:
	@ if [ ! -d "${${*}}" ]; then \
		echo "Directory $* does not exist"; \
		exit 1; \
	fi


TOPHAT_OUTPUT_REL ?= "tophat_output"
CUFFLINKS_OUTPUT_REL ?= "cufflinks_output"

cufflinks-serial: guard-env-LIBRARY_DIR guard-dir-LIBRARY_DIR
	@TOPHAT_OUTPUT_DIR="${LIBRARY_DIR}/${TOPHAT_OUTPUT_REL}" \
	CUFFLINKS_OUTPUT_DIR="${LIBRARY_DIR}/${CUFFLINKS_OUTPUT_REL}" \
	./cufflinks-serial.sh


queue-jobs: guard-env-LIBRARY_DIR guard-dir-LIBRARY_DIR
	@TOPHAT_OUTPUT_DIR="${LIBRARY_DIR}/${TOPHAT_OUTPUT_REL}" \
	CUFFLINKS_OUTPUT_DIR="${LIBRARY_DIR}/${CUFFLINKS_OUTPUT_REL}" \
	./queue-jobs.sh
