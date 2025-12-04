# Burke Lab Genomics Pipeline
## Powered by Nextflow

This is a README file for the Burke Lab Genomics Pipeline. It's purpose is to explain what's going on with the pipeline and how to run it properly.

The purpose of this pipeline is to create a SNP and INDEL table from genomic data of many samples.
You can start with a number of .fastq files for a number of samples OR one intermediate vcf file per sample.

This pipeline follows this organizational structure:
```
ProjectName (e.g. 2023_FLYLONG)/
├── raw
│   ├── sample1_R1.fastq.gz
│   ├── sample1_R2.fastq.gz
│   ├── sample2_R1.fastq.gz
│   └── [...]
├── vcf_files/
│   ├── pop1/
│   │   ├── pop1.vcf
│   │   └── pop1.vcf.tbi
│   └── [...]
└── pipeline_runs/
    ├── 2025_01_01_round01
    │   ├── results/
    │   │   ├── filtered_snps.txt
    │   │   └── [...]
    │   ├── logs/
    │   │   └── [...]
    │   ├── reports/
    │   │   ├── fastQC/
    │   │   │   ├── pop1/
    │   │   │   │   └── [...] 
    │   │   │   └── [...]
    │   │   ├── coverage/
    │   │   │   ├── pop1/
    │   │   │   │   └── [...] 
    │   │   │   └── [...]
    │   │   └── [...]
    │   └── metadata.csv
    └── 2025_02_01_round02
```

## How to run the pipeline
1. First, ensure that you have your project directory (inside the lab's folder) with a subdirectory `raw` containing one subdirectory for each of your samples, and in each sample subdirectory, their respective `.fastq` files. Alternatively, you can start from intermediate vcf files. To do that, ensure that your project directory have a `vcf_files` subdirectory containing one subdirectory per sample, and in each subdirectory the `.vcf` (and `.tbi`) file of that sample.
2. Inside your project directory, create the subdirectories `pipeline_runs`. Inside `pipeline_runs`, create a subdirectory with the name of your pipeline run. I advise it to be 'YYYY_MM_DD_roundXX'. That keeps the subdirectories organized and listed chronologically.
3. Make sure that you have the `metadata_nextflow.csv` file inside the pipeline run directory.
4. When that is ready, you can run the pipeline. Edit `nextflow.config`, changing `params.project` and `params.pipeline_run` with the names of your project (e.g. `2025_MULTI`) and pipeline run name (e.g. `2025_01_01_round01`).
5. Cleanup time - run these steps to ensure that your scratch space is clean:
```
salloc -p burke_lab
rm -rf /scratch/nextflow_pipeline_scratch/*
exit
```
6. Next, make sure `run.sh` is executable (`chmod +x run.sh`), open it and edit `RUN_ID` to be `attempt_01`. Sometimes the pipeline run doesn't go as we expect. So every time you have to re-run the pipeline, you change `RUN_ID` (`attempt_02`, `attempt_03`, etc). That will work only if the nextflow command has the `-resume` flag. You may want to remove it if you want to run things anew.
7. You can run `hqstat --watch` to see your pipeline run's progress. If it says `running`, its working! If it says `failed`, check your pipeline run logs (it's in the directory just created, `attempt_XX`, usually `attempt_XX.oNNNNNNN`), fix the problems, change `RUN_ID` inside `run.sh`, and try again!