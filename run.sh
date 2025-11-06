#!/bin/bash

echo "This will submit all the jobs. Submit? (yes/no): "
read input
if [ "$input" != "yes" ]; then
  echo "Quitting."
  exit
fi

# Run pipeline
echo "nextflow run main.nf -ansi-log false" | hqsub - -q burke_lab -t array -r attempt_24 --local-drive shared --local-prefix /scratch/workingdir


