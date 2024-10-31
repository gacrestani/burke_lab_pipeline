#!/bin/bash

# Run pipeline
echo "nextflow run main.nf -ansi-log false" | hqsub - -q burke_lab -t array -r gatk_pipeline_01 --local-drive shared --local-prefix /scratch
