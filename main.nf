log.info """\
    B U R K E   L A B   P I P E L I N E
    ===================================
    """
    .stripIndent()

process BwaMem {
    label 'low_reqs'
    debug params.debug
    
    input:
    tuple val(population), val(meta), path(reads)

    output:
    tuple val(population), path("${reads[0].simpleName}_sorted.bam"), emit: bam

    script:
    def fastqs = reads.collect().join(' ')
    """
    # bwa mem script
    # Defining the bwa-mem/samtools command

    # Defining RG
    RG="@RG\\tID:${meta.flow_cell}.lane-${meta.lane}.${meta.barcode}\\tSM:${population}\\tLB:${meta.internal_library_name}\\tPL:ILLUMINA\\tPU:${meta.flow_cell}.${meta.lane}.${meta.barcode}"

    # Defining the bwa mem | samtools command
    cmd="bwa mem -R '\$RG' ${params.reference_genome} -t ${task.cpus} ${fastqs} | samtools sort --threads ${task.cpus} -o ${reads[0].simpleName}_sorted.bam"
    
    # Logging command
    echo "\$cmd"

    eval \$cmd
    """
    
    stub:
    """
    touch ${reads[0].simpleName}_sorted.bam
    """
}

process BamToVcf {
    label 'medium_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.txt"
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.table"
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.bai"
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.g.vcf.gz"

    input:
    tuple val(population), path(bams), val(bams_path)

    output:
    val(population)
    path("${population}_haplotype_caller_results.g.vcf.gz"), emit: vcf
    path("${population}_haplotype_caller_results.g.vcf.gz.tbi"), emit: vcftbi

    script:
    def bams_list = bams.collect{"--INPUT $it"}.join(' ')
    """
    echo "Starting BamToVcf..."
    echo "Space used before the script:"
    df -h .
    echo "Space used in the scratch directory:"
    df -h ${params.scratch_directory}

    ############################
    # gatk MergeSamFiles script
    # Defining the command
    cmd="gatk MergeSamFiles ${bams_list} --OUTPUT mergesamfiles.bam --TMP_DIR ${params.scratch_directory}"

    echo "\$cmd"

    # Run command
    eval \$cmd
    ############################


    ############################
    # gatk MarkDuplicates script
    # Defining the command
    cmdMD="gatk MarkDuplicates --INPUT mergesamfiles.bam --METRICS_FILE duplicate_metrics.txt --OUTPUT duplicates_marked.bam --TMP_DIR ${params.scratch_directory}"

    echo "\$cmdMD"
   
    # Run command
    eval \$cmdMD
    ############################

    ############################
    # gatk BaseRecalibrator script
    # Defining the command
    cmdBR="gatk BaseRecalibrator --input duplicates_marked.bam --known-sites ${params.bqsr_vcf} --output recalibration_metrics.table --reference ${params.reference_genome} --tmp-dir ${params.scratch_directory}"

    echo "\$cmdBR"
    
    # Run command
    eval \$cmdBR
    ############################


    ###########################
    # gatk ApplyBQSR script
    # Defining the command
    cmdABQSR="gatk ApplyBQSR --bqsr-recal-file recalibration_metrics.table --input duplicates_marked.bam --output haplotype_this.bam --tmp-dir ${params.scratch_directory}"

    echo "\$cmdABQSR"
    
    # Run command
    eval \$cmdABQSR
    ###########################


    ##########################
    # gatk HaplotypeCaller script
    # Defining the command
    cmdHC="gatk HaplotypeCaller -R ${params.reference_genome} -ERC GVCF -I haplotype_this.bam -O ${population}_haplotype_caller_results.g.vcf.gz --tmp-dir ${params.scratch_directory}"

    echo "\$cmdHC"
    
    # Run command
    eval \$cmdHC
    ###########################

    echo "Cleanup all .bam files"
    rm -rf *.bam

    echo "Space used after script is done:"
    du -sh .

    echo "Space used in the scratch directory:"
    du -sh ${params.scratch_directory}
    """

    stub:
    """
    touch duplicate_metrics.txt
    touch recalibration_metrics.table
    touch file.bai

    touch ${population}_haplotype_caller_results.g.vcf.gz
    touch ${population}_haplotype_caller_results.g.vcf.gz.tbi
    """
}


process MergeSamFiles {
    label 'low_reqs'
    debug params.debug

    input:
    tuple val(population), path(bams)

    output:
    tuple val(population), path("${population}.bam"), emit: bam

    script:
    def bams_list = bams.collect{"--INPUT $it"}.join(' ')
    """
    # gatk MergeSamFiles script
    # Defining the command
    cmd="gatk MergeSamFiles ${bams_list} --OUTPUT ${population}.bam --TMP_DIR ${params.scratch_directory}"

    echo "\$cmd"

    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${population}.bam
    """
}

process MarkDuplicates {
    label 'medium_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(population), path(bam)

    output:
    tuple val(population), path("${population}_duplicates_marked.bam"), emit: bam
    path("${population}_duplicate_metrics.txt")

    script:
    """
    # gatk MarkDuplicates script
    # Defining the command
    cmd="gatk MarkDuplicates --INPUT ${bam} --METRICS_FILE ${population}_duplicate_metrics.txt --OUTPUT ${population}_duplicates_marked.bam --TMP_DIR ${params.scratch_directory}"

    echo "\$cmd"
   
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${population}_duplicates_marked.bam
    touch ${population}_duplicate_metrics.txt
    """
}

process BaseRecalibrator {
    label 'medium_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.table"

    input:
    tuple val(population), path(duplicates_marked_bam)

    output:
    tuple val(population), path(duplicates_marked_bam), path("${population}_recalibration_metrics.table"), emit: bam

    script:
    """
    # gatk BaseRecalibrator script
    # Defining the command
    cmd="gatk BaseRecalibrator --input ${duplicates_marked_bam} --known-sites ${params.bqsr_vcf} --output ${population}_recalibration_metrics.table --reference ${params.reference_genome} --tmp-dir ${params.scratch_directory}"

    echo "\$cmd"
    
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${population}_duplicates_marked.bam
    touch ${population}_recalibration_metrics.table
    """
}

process ApplyBQSR {
    label 'medium_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.bai"

    input:
    tuple val(population), path(duplicates_marked_bam), path(recalibration_metrics_table)

    output:
    tuple val(population), path("${population}_haplotype_this.bam"), path("*.bai"), emit: bam

    script:
    """
    # gatk ApplyBQSR script
    # Defining the command
    cmd="gatk ApplyBQSR --bqsr-recal-file ${recalibration_metrics_table} --input ${duplicates_marked_bam} --output ${population}_haplotype_this.bam --tmp-dir ${params.scratch_directory}"

    echo "\$cmd"
    
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${population}_haplotype_this.bam
    touch ${population}_haplotype_this.bam.bai
    """
}

process HaplotypeCaller {
    label 'medium_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/${population}", mode: 'copy', pattern: "*.g.vcf.gz"

    input:
    tuple val(population), path(haplotype_this_bam), path(bai)

    output:
    val(population)
    path("${population}_haplotype_caller_results.g.vcf.gz"), emit: vcf
    path("${population}_haplotype_caller_results.g.vcf.gz.tbi"), emit: vcftbi

    script:
    """
    # gatk HaplotypeCaller script
    # Defining the command
    cmd="gatk HaplotypeCaller -R ${params.reference_genome} -ERC GVCF -I ${haplotype_this_bam} -O ${population}_haplotype_caller_results.g.vcf.gz --tmp-dir ${params.scratch_directory}"

    echo "\$cmd"
    
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${population}_haplotype_caller_results.g.vcf.gz
    touch ${population}_haplotype_caller_results.g.vcf.gz.tbi
    """
}

process CombineGVCFs {
    label 'high_reqs'
    debug params.debug

    input:
    path(haplotype_caller_results)
    path(haplotype_caller_results_tbi)

    output:
    tuple path("combined.g.vcf.gz"), path("combined.g.vcf.gz.tbi"), emit: vcf

    script:
    def variant_list = haplotype_caller_results.collect{"--variant $it"}.join(' ')
    def variant_tbi_list = haplotype_caller_results_tbi.collect{"--read-index $it"}.join(' ')

    """
    # gatk CombineGVCFs script
    # Defining the command
    cmd="gatk CombineGVCFs -R ${params.reference_genome} ${variant_list} -O combined.g.vcf.gz ${variant_tbi_list} --tmp-dir ${params.scratch_directory}"

    echo "\$cmd"
    
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch combined.g.vcf.gz
    touch combined.g.vcf.gz.tbi
    """
}

process GenotypeGVCFs {
    label 'high_reqs'
    debug params.debug

    input:
    tuple path(combined_g_vcf_gz), path(combined_g_vcf_gz_tbi)

    output:
    tuple path("all_samples_raw.vcf.gz"), path("all_samples_raw.vcf.gz.tbi"), emit: vcf

    script:
    """
    # gatk GenotypeGVCFs script
    # Defining the command
    cmd="gatk GenotypeGVCFs -R ${params.reference_genome} -V ${combined_g_vcf_gz} -O all_samples_raw.vcf.gz --read-index ${combined_g_vcf_gz_tbi} --tmp-dir ${params.scratch_directory}"

    echo "\$cmd"
    
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch all_samples_raw.vcf.gz
    touch all_samples_raw.vcf.gz.tbi
    """
}

process SelectVariants {
    label 'high_reqs'
    debug params.debug

    input:
    tuple path(all_samples_raw_vcf_gz), path(all_samples_raw_vcf_gz_tbi)

    output:
    tuple path("raw_snps.vcf.gz"), path("raw_indels.vcf.gz"), path("raw_snps.vcf.gz.tbi"), path("raw_indels.vcf.gz.tbi"), emit: vcf

    script:
    """
    # gatk SelectVariants script
    # Defining commands
    cmdA="gatk SelectVariants --variant ${all_samples_raw_vcf_gz} --output raw_snps.vcf.gz --select-type-to-include SNP --read-index ${all_samples_raw_vcf_gz_tbi} --tmp-dir ${params.scratch_directory}"
    cmdB="gatk SelectVariants --variant ${all_samples_raw_vcf_gz} --output raw_indels.vcf.gz --select-type-to-include INDEL --select-type-to-include MIXED --exclude-non-variants --read-index ${all_samples_raw_vcf_gz_tbi} --tmp-dir ${params.scratch_directory}"
    
    echo "\$cmdA"
    echo "\$cmdB"

    # Run command
    eval \$cmdA
    eval \$cmdB
    """

    stub:
    """
    touch raw_snps.vcf.gz
    touch raw_snps.vcf.gz.tbi
    touch raw_indels.vcf.gz
    touch raw_indels.vcf.gz.tbi
    """
}

process VariantFiltration {
    label 'high_reqs'
    debug params.debug

    input:
    tuple path(raw_snps_vcf_gz), path(raw_indels_vcf_gz), path(raw_snps_vcf_gz_tbi), path(raw_indels_vcf_gz_tbi)

    output:
    tuple path("filtered_snps.vcf"), path("filtered_indels.vcf"), emit: vcf

    script:
    """
    # gatk VariantFiltration script
    # Defining command
    cmdA="gatk VariantFiltration --read-index ${raw_snps_vcf_gz_tbi} --variant ${raw_snps_vcf_gz} --output filtered_snps.vcf --filter-name LowVQCBD --filter-expression \\"QD < 5.0\\" --filter-name LowQualSNPs --filter-expression \\"MQ < 40.0\\"  --filter-name FisherStrand --filter-expression \\"FS > 60.0\\" --filter-name ReadPosRank --filter-expression \\"ReadPosRankSum < -8.0\\" --filter-name MQRankSum --filter-expression \\"MQRankSum < -12.5\\" --tmp-dir ${params.scratch_directory}"
    cmdB="gatk VariantFiltration --read-index ${raw_indels_vcf_gz_tbi} --variant ${raw_indels_vcf_gz} --output filtered_indels.vcf  --filter-name LowVQCBD --filter-expression \\"QD < 5.0\\" --filter-name FisherStrand --filter-expression \\"FS > 200.0\\" --filter-name ReadPosRank --filter-expression \\"ReadPosRankSum < -20.0\\" --tmp-dir ${params.scratch_directory}"

    echo "\$cmdA"
    echo "\$cmdB"

    # Run command
    eval \$cmdA
    eval \$cmdB
    """

    stub:
    """
    touch filtered_snps.vcf
    touch filtered_indels.vcf
    """
}

process SnpEff {
    label 'high_reqs'
    debug params.debug

    input:
    path(filtered_vcf)

    output:
    path("${filtered_vcf.baseName}_ann.vcf"), emit: vcf

    script:
    """
    # SnpEff script
    # Defining commands
    cmd="snpEff -v -o gatk -dataDir ${params.scratch_directory} ${params.snpeff_organism} ${filtered_vcf} > ${filtered_vcf.baseName}_ann.vcf"

    echo "\$cmd"

    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${filtered_vcf.baseName}_ann.vcf
    """
}

process VariantsToTable {
    label 'high_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/", mode: 'copy', pattern: "filtered_snps_ann.txt", saveAs: { filename -> "annotated_snps.txt" }
    publishDir path: "${params.results_directory}/", mode: 'copy', pattern: "filtered_indels_ann.txt", saveAs: { filename -> "annotated_indels.txt" }

    input:
    path(filtered_ann_vcf)

    output:
    path("${filtered_ann_vcf.baseName}.txt"), emit: txt

    script:
    """
    # gatk VariantsToTable script
    # Defining command
    cmd="gatk VariantsToTable -V ${filtered_ann_vcf} -F CHROM -F POS -F REF -F ALT -F HOM-VAR -F HET -F HOM-REF -F NO-CALL -F TYPE -F EVENTLENGTH -F MULTI-ALLELIC -F ANN -GF GT -GF AD -O ${filtered_ann_vcf.baseName}.txt --tmp-dir ${params.scratch_directory}"

    echo "\$cmd"

    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${filtered_ann_vcf.baseName}.txt
    touch filtered_snps_ann.txt
    touch filtered_indels_ann.txt
    """
}

process VcfToTable {
    label 'high_reqs'
    debug params.debug
    publishDir path: "${params.results_directory}/", mode: 'copy', pattern: "*.txt"

    input:
    path(filtered_vcf)

    output:
    path("${filtered_vcf.baseName}.txt"), emit: txt

    script:
    """
    # VcfToTable script
    # Defining the command
    cmd="/nfs3/IB/Burke_Lab/Crestani/nextflow/local/bin/vcf_to_table.py --vcf ${filtered_vcf} --output ${filtered_vcf.baseName}.txt --num_allow_missing 0"

    echo "\$cmd"
    
    # Run command
    eval \$cmd
    """

    stub:
    """
    touch ${filtered_vcf.baseName}.txt
    """
}


workflow {
    Channel.fromPath(params.metadata)
    | splitCsv(header:true)
    | map { row ->
        fastq1_path = params.samples_directory + row.fastq1
        fastq2_path = row.fastq2 ? params.samples_directory + row.fastq2 : null

        meta = row.subMap(
            'flow_cell',
            'lane',
            'population',
            'barcode',
            'sequencing_facility',
            'internal_library_name'
            )

        def files = [file(fastq1_path, checkIfExists: true)]
        if (fastq2_path) {
            files << file(fastq2_path, checkIfExists: true)
        }

        [row.population, meta, files]
    }
//    | filter { it.contains("EB_rep04_gen20") || it.contains("CB_rep01_gen56")}
    | set { samples }

    BwaMem(samples)

    // Checks if params.run_serialized is TRUE. If so, run the process that takes less scratch space.
    if( params.run_serialized ) {

        BamToVcf(BwaMem.out.bam.groupTuple())
        vcf_ch = BamToVcf.out

    } else {

        MergeSamFiles(BwaMem.out.bam.groupTuple())
        MarkDuplicates(MergeSamFiles.out.bam)
        BaseRecalibrator(MarkDuplicates.out.bam)
        ApplyBQSR(BaseRecalibrator.out.bam)
        HaplotypeCaller(ApplyBQSR.out.bam)
        vcf_ch = HaplotypeCaller.out

    }

    CombineGVCFs(vcf_ch.vcf.collect(), vcf_ch.vcftbi.collect())
    GenotypeGVCFs(CombineGVCFs.out.vcf)
    SelectVariants(GenotypeGVCFs.out.vcf)
    VariantFiltration(SelectVariants.out.vcf)
    
    SnpEff(VariantFiltration.out.vcf.flatten())
    VariantsToTable(SnpEff.out.vcf)
    
    VcfToTable(VariantFiltration.out.vcf.flatten())
}
