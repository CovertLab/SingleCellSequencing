.PHONY: cufflinks-serial queue-jobs

guard-%:
	@ if [ "${${*}}" == "" ]; then \
		echo "Environment variable $* not set"; \
		exit 1; \
	fi

TOPHAT_OUTPUT_REL ?= "tophat_output"
CUFFLINKS_OUTPUT_REL ?= "cufflinks_output"

cufflinks-serial: guard-LIBRARY_DIR
	@TOPHAT_OUTPUT_DIR="${LIBRARY_DIR}/${TOPHAT_OUTPUT_REL}" \
	CUFFLINKS_OUTPUT_DIR="${LIBRARY_DIR}/${CUFFLINKS_OUTPUT_REL}" \
	./cufflinks-serial.sh


queue-jobs: guard-LIBRARY_DIR
	@TOPHAT_OUTPUT_DIR="${LIBRARY_DIR}/${TOPHAT_OUTPUT_REL}" \
	CUFFLINKS_OUTPUT_DIR="${LIBRARY_DIR}/${CUFFLINKS_OUTPUT_REL}" \
	./queue-jobs.sh
