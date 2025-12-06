#!/bin/bash

RUN_ID="attempt_04"








echo "This will submit all the jobs. Submit? (yes|y / no): "
read input
if [[ "$input" != "yes" && "$input" != "y" ]]; then
  echo "Quitting."
  exit
fi

# Run pipeline
echo "nextflow run main.nf -ansi-log false -resume" | hqsub - -q burke_lab -t array -r "$RUN_ID" --local-drive shared --local-prefix /scratch/nextflow_pipeline_scratch


