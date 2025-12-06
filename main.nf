log.info """\
    B U R K E   L A B   G E N O M I C   D A T A   P I P E L I N E
    =============================================================
    This is the Burke Lab Genomic Pipeline main script.
    There's no need to edit anything here! Go to nextflow.config
    to set the parameters for your pipeline run. Good luck!
    """
    .stripIndent()

// Processes ===================================================================

process BwaMem {
    // This process performs BWA MEM alignment on the provided read files
    // It then sorts the resulting SAM file into a BAM file
    // It takes the reads and meta information, runs the alignment, and outputs a sorted BAM file

    label 'parallel_per_sample' 
    debug params.debug
    
    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${population}/${task.hash}",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: population, metadata, and read files
    input: 
    tuple val(population), val(meta), path(reads), path(vcf) 
    
    // Output: Sorted BAM file
    output: 
    tuple val(population), path("${reads[0].simpleName}_sorted.bam"), emit: bam
    path(".command.*"), emit: logs

    script:
    // Define the Read Group string
    def rg_string = "@RG\\tID:${meta.flow_cell}.${meta.lane}\\tSM:${population}\\tLB:${meta.internal_library_name}\\tPL:ILLUMINA\\tPU:${meta.flow_cell}.${meta.lane}.${meta.barcode}"

    // Efficiently split allocated CPUs between the two piped commands
    def bwa_threads = (task.cpus * 0.75).toInteger() ?: 1
    def samtools_threads = (task.cpus * 0.25).toInteger() ?: 1
    
    // Set memory per thread for samtools sort, ensuring a minimum of 2G
    def mem_val = (task.memory.toGiga() / samtools_threads).toInteger()
    def mem_per_thread = "${Math.max(mem_val, 2)}G"
    
    """
    bwa mem -M \\
        -t ${bwa_threads} \\
        -R '${rg_string}' \\
        ${params.reference_genome} \\
        ${reads.join(' ')} \\
    | samtools sort \\
        --threads ${samtools_threads} \
        -m ${mem_per_thread} \\
        -o ${reads[0].simpleName}_sorted.bam
    """
    
    stub:
    """
    touch ${reads[0].simpleName}_sorted.bam
    """
}

process MergeSamFiles {
    // This process merges multiple BAM files into a single BAM file using GATK's MergeSamFiles tool
    // It takes several BAM files as input, and outputs the merged BAM file.

    label 'io_intensive'
    debug params.debug

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${population}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: population and list of BAM files
    input:
    tuple val(population), path(bam)

    // Output: Merged BAM file, and its index BAI file
    output:
    tuple val(population), path("${population}.bam"), path("${population}.bai"), emit: bam
    path(".command.*"), emit: logs

    // First create a bams_list with all bam files passed as input
    script:
    def bams_list = bam.collect{"--INPUT $it"}.join(' ')
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \\
        MergeSamFiles \\
        ${bams_list} \\
        --OUTPUT ${population}.bam \\
        --CREATE_INDEX true \\
        --USE_THREADING true \\
        --TMP_DIR ${params.scratch_directory}
    """

    stub:
    """
    touch ${population}.bam
    touch ${population}.bai
    """
}

process MarkDuplicates {
    // This process marks duplicates in a BAM file using Picard's MarkDuplicates tool
    // It takes a BAM file as input, marks the duplicates, and outputs a BAM file with duplicates marked.

    label 'high_mem_single'
    debug params.debug

    // Publishing the duplicate_metrics.txt file if required for future analyses
    publishDir (
        path: "${params.results_directory}/intermediate_files/${population}/",
        mode: 'copy',
        pattern: "*.txt"
    )
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${population}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: population, BAM file, and its index BAI file
    input:
    tuple val(population), path(bam), path(bai)

    // Output: BAM file with duplicates marked, and its index BAI file
    output:
    tuple (val(population), path("${population}_duplicates_marked.bam"), path("${population}_duplicates_marked.bai"), emit: bam)
    path("${population}_duplicate_metrics.txt")
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \\
        MarkDuplicates \\
        --INPUT ${bam} \\
        --METRICS_FILE ${population}_duplicate_metrics.txt \\
        --OUTPUT ${population}_duplicates_marked.bam \\
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
        --VALIDATION_STRINGENCY SILENT \\
        --ASSUME_SORTED true \\
        --CREATE_INDEX true \\
        --TMP_DIR ${params.scratch_directory}
    """

    stub:
    """
    touch ${population}_duplicates_marked.bam
    touch ${population}_duplicates_marked.bai
    touch ${population}_duplicate_metrics.txt
    """
}

process BaseRecalibrator {
    // This process runs GATK's BaseRecalibrator to recalibrate base quality scores
    // It takes a sorted BAM file and a reference genome, applies base recalibration,
    // and produces a recalibrated BAM file and a recalibration report.

    label 'high_mem_single'
    debug params.debug

    // Publishing the recalibration table if required for future analyses
    publishDir (
        path: "${params.results_directory}/intermediate_files/${population}",
        mode: 'copy',
        pattern: "*.table"
    )
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${population}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: population, duplication-marked BAM file, and its index TBI file
    // This process uses the known variant sites (reference VCF file)
    input:
    tuple val(population), path(duplicates_marked_bam), path(duplicates_marked_bai)

    // Output: recalibrated BAM and recalibration report
    output:
    tuple val(population), path(duplicates_marked_bam), path(duplicates_marked_bai), path("${population}_recalibration_metrics.table"), emit: bam
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \\
        BaseRecalibrator \\
        --input ${duplicates_marked_bam} \\
        --known-sites ${params.bqsr_vcf} \\
        --output ${population}_recalibration_metrics.table \\
        --reference ${params.reference_genome} \\
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch ${population}_recalibration_metrics.table
    """
}

process ApplyBQSR {
    // This process applies base quality score recalibration to a BAM file using a recalibration report
    // It takes the original BAM file and a recalibration report, applies the recalibration,
    // and outputs the recalibrated BAM file.

    label 'high_mem_single'
    debug params.debug
    
    // Publishing the final BAM Files before turning then into VCFs
    // Change: no need to save .bam files. Commenting out this section.
    // publishDir (
    //     path: "${params.results_directory}/intermediate_files/${population}",
    //     mode: 'copy',
    //     pattern: "*.bam*" // Publishes both .bam and .bai
    // )
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${population}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: population, duplication-marked BAM file, reference genome, and recalibration report
    input:
    tuple val(population), path(duplicates_marked_bam), path(duplicates_marked_bai), path(recalibration_metrics_table)

    // Output: Recalibrated BAM file
    output:
    tuple val(population), path("${population}_haplotype_this.bam"), path("${population}_haplotype_this.bai"), emit: bam
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \\
        ApplyBQSR \\
        --reference ${params.reference_genome} \\
        --bqsr-recal-file ${recalibration_metrics_table} \\
        --create-output-bam-index true \\
        --input ${duplicates_marked_bam} \\
        --output ${population}_haplotype_this.bam \\
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch ${population}_haplotype_this.bam
    touch ${population}_haplotype_this.bai
    """
}

process HaplotypeCaller {
    // This process runs GATK's HaplotypeCaller to call variants from a recalibrated BAM file
    // It takes a recalibrated BAM file and a reference genome, applies the variant calling,
    // and outputs a GVCF file with the called variants.

    label 'parallel_per_sample'
    debug params.debug

    // Publishing the VCF file and its index TBI file
    publishDir (
        path: "${params.vcf_directory}/${population}",
        mode: 'copy',
        pattern: "*.g.vcf.gz*" // Catches both .gz and .tbi
    )
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${population}/",
        mode: 'copy',
        pattern: ".command.*"
    )
    
    // Input: population, recalibrated BAM file, reference genome
    input:
    tuple val(population), path(haplotype_this_bam), path(haplotype_this_bai)

    // Output: G.VCF file and its index TBI file
    output:
    path("${population}_haplotype_caller_results.g.vcf.gz"), emit: vcf
    path("${population}_haplotype_caller_results.g.vcf.gz.tbi"), emit: tbi
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \\
        HaplotypeCaller \\
        -R ${params.reference_genome} \\
        -I ${haplotype_this_bam} \\
        -O ${population}_haplotype_caller_results.g.vcf.gz \\
        -ERC GVCF \\
        --native-pair-hmm-threads ${task.cpus} \\
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch ${population}_haplotype_caller_results.g.vcf.gz
    touch ${population}_haplotype_caller_results.g.vcf.gz.tbi
    """
}

process CombineGVCFs { // Maybe GenomicsDB is better here
    // This process combines multiple GVCF files into a single GVCF file using GATK's CombineGVCFs tool
    // It takes multiple GVCF files as input, combines them, and outputs a single merged GVCF file.

    label 'joint_calling'
    debug params.debug

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: list of GVCF files and their index TBI files
    // Control here serves as a way to stop this from happening before all HaplotypeCaller calls end
    input:
    path(vcfs)
    path(tbis)

    // Output: Combined GVCF file and its index TBI file
    output:
    tuple (
        path("combined.g.vcf.gz"),
        path("combined.g.vcf.gz.tbi"),
        emit: vcf
    )
    path(".command.*"), emit: logs

    script:
    def variant_list = vcfs.collect{"--variant ${it}"}.join(' ')
    def variant_tbi_list = tbis.collect{"--read-index ${it}"}.join(' ')
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \
        CombineGVCFs \
        -R ${params.reference_genome} \
        ${variant_list} \
        -O combined.g.vcf.gz \
        ${variant_tbi_list} \
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch combined.g.vcf.gz
    touch combined.g.vcf.gz.tbi
    """
}

process GenotypeGVCFs {
    // This process runs GATK's GenotypeGVCFs to perform variant calling on a combined GVCF file
    // It takes a combined GVCF file, applies the genotyping, and outputs a VCF file with the called variants.

    label 'joint_calling'
    debug params.debug

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: combined GVCF file and its index TBI file
    input:
    tuple (
        path(combined_g_vcf_gz),
        path(combined_g_vcf_gz_tbi)
    )

    // Output: VCF file containing the variant calls and its index TBI file
    output:
    tuple (
        path("all_samples_raw.vcf.gz"),
        path("all_samples_raw.vcf.gz.tbi"),
        emit: vcf
    )
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \
        GenotypeGVCFs \
        -R ${params.reference_genome} \
        -V ${combined_g_vcf_gz} \
        -O all_samples_raw.vcf.gz \
        --read-index ${combined_g_vcf_gz_tbi} \
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch all_samples_raw.vcf.gz
    touch all_samples_raw.vcf.gz.tbi
    """
}

process SelectVariants { // This could be broken into two different processes
    // This process runs GATK's SelectVariants to filter and select specific variants from a VCF file
    // It takes a VCF file, applies selection criteria (e.g., specific variants, regions, etc.),
    // and outputs a filtered VCF file containing only the selected variants.

    label 'high_mem_single'
    debug params.debug

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/",
        mode: 'copy', pattern: ".command.*"
    )

    // Input: genotyped VCF file and its index TBI file
    input:
    tuple (
        path(all_samples_raw_vcf_gz),
        path(all_samples_raw_vcf_gz_tbi)
    )

    // Output: Filtered VCF file containing the selected variants, and their index TBI files
    output:
    tuple (
        path("raw_snps.vcf.gz"),
        path("raw_indels.vcf.gz"),
        path("raw_snps.vcf.gz.tbi"),
        path("raw_indels.vcf.gz.tbi"),
        emit: vcf
    )
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \
        SelectVariants \
        --variant ${all_samples_raw_vcf_gz} \
        --output raw_snps.vcf.gz \
        --select-type-to-include SNP \
        --read-index ${all_samples_raw_vcf_gz_tbi} \
        --tmp-dir ${params.scratch_directory}

    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \
        SelectVariants \
        --variant ${all_samples_raw_vcf_gz} \
        --output raw_indels.vcf.gz \
        --select-type-to-include INDEL \
        --select-type-to-include MIXED \
        --exclude-non-variants \
        --read-index ${all_samples_raw_vcf_gz_tbi} \
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch raw_snps.vcf.gz
    touch raw_snps.vcf.gz.tbi
    touch raw_indels.vcf.gz
    touch raw_indels.vcf.gz.tbi
    """
}

process VariantFiltration { // Could be broken into two processes
    // This process runs GATK's VariantFiltration to apply filtering criteria to a VCF file
    // It takes a VCF file, applies the specified filters (e.g., quality score, depth),
    // and outputs a filtered VCF file with only the variants passing the criteria.

    label 'high_mem_single'
    debug params.debug

    // Publishing the filtered VCF
    publishDir (
        path: "${params.results_directory}",
        mode: 'copy', pattern: "*.vcf"
    )

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: VCF files and their index TBI files
    input:
    tuple (
        path(raw_snps_vcf_gz),
        path(raw_indels_vcf_gz),
        path(raw_snps_vcf_gz_tbi),
        path(raw_indels_vcf_gz_tbi)
    )

    // Output: Filtered VCF file containing only variants passing the filtering criteria, one for SNPs and one for Indels
    output:
    tuple (
        path("filtered_snps.vcf"),
        path("filtered_indels.vcf"),
        emit: vcf
    )
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \
        VariantFiltration \
        --read-index ${raw_snps_vcf_gz_tbi} \
        --variant ${raw_snps_vcf_gz} \
        --output filtered_snps.vcf \
        --filter-name LowVQCBD \
        --filter-expression 'QD < 5.0' \
        --filter-name LowQualSNPs \
        --filter-expression 'MQ < 40.0'  \
        --filter-name FisherStrand \
        --filter-expression 'FS > 60.0' \
        --filter-name ReadPosRank \
        --filter-expression 'ReadPosRankSum < -8.0' \
        --filter-name MQRankSum \
        --filter-expression 'MQRankSum < -12.5' \
        --tmp-dir ${params.scratch_directory}

    ${params.gatk} --java-options '-Xmx${task.memory.giga}G' \
        VariantFiltration \
        --read-index ${raw_indels_vcf_gz_tbi} \
        --variant ${raw_indels_vcf_gz} \
        --output filtered_indels.vcf  \
        --filter-name LowVQCBD \
        --filter-expression 'QD < 5.0' \
        --filter-name FisherStrand \
        --filter-expression 'FS > 200.0' \
        --filter-name ReadPosRank \
        --filter-expression 'ReadPosRankSum < -20.0' \
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch filtered_snps.vcf
    touch filtered_indels.vcf
    """
}

process SnpEff {
    // This process runs SnpEff to annotate variants in a VCF file
    // It takes a VCF file, annotates the variants using a provided genome database,
    // and outputs a VCF file with the annotations added.

    // WARNING: snpEff is not working well with yeast data.
    // I believe it has to do with the chromosome renaming we have to do
    // we use chromosome1 instead of chrom1 (or vice-versa)

    label 'high_mem_single'
    debug params.debug


    // Publishing all files generated by SnpEff
    publishDir (
        path: "${params.results_directory}/",
        mode: 'copy',
        pattern: '*_ann.vcf'
    )
    publishDir (
        path: "${params.results_directory}/reports/${task.process}/${task.hash}",
        mode: 'copy',
        pattern: '*.txt'
    )
    publishDir (
        path: "${params.results_directory}/reports/${task.process}/${task.hash}",
        mode: 'copy',
        pattern: '*.html'
    )
    publishDir (
        path: "${params.results_directory}/reports/${task.process}/${task.hash}",
        mode: 'copy',
        pattern: '*.csv'
    )
    
    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/${task.hash}",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: VCF file
    input:
    path(filtered_vcf)

    // Output: annotated VCF file
    output:
    path("${filtered_vcf.baseName}_ann.vcf"), emit: vcf
    path(".command.*"), emit: logs

    script:
    def heap_mem = (task.memory.toGiga() * 0.9 * 0.5).toInteger()
    """
    java -Xmx${heap_mem}G -Xms${heap_mem}G -jar /nfs6/TMP/Burke_Lab/shared_data/snpEff/snpEff.jar \
        -v -dataDir ${params.scratch_directory} \
        -csvStats snpeff_stats.csv \
        ${params.snpeff_organism} \
        ${filtered_vcf} \
        > ${filtered_vcf.baseName}_ann.vcf
    """

    stub:
    """
    touch ${filtered_vcf.baseName}_ann.vcf
    """
}

process VariantsToTable {
    // This process runs GATK's VariantsToTable to convert a VCF file into a table format
    // Currently being used to create tables from the files annotated by SnpEff
    // This output file could work as our primary snp_table too
    // It extracts specified fields from the VCF file and outputs them in a tabular format.

    label 'io_intensive'
    debug params.debug

    // Publishing the txt files generated
    publishDir (
        path: "${params.results_directory}/",
        mode: 'copy',
        pattern: "filtered_snps_ann.txt",
        saveAs: { filename -> "annotated_snps.txt" }
    )
    publishDir (
        path: "${params.results_directory}/",
        mode: 'copy',
        pattern: "filtered_indels_ann.txt",
        saveAs: { filename -> "annotated_indels.txt" }
    )

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    input:
    path(filtered_ann_vcf)

    output:
    path("${filtered_ann_vcf.baseName}.txt"), emit: txt
    path(".command.*"), emit: logs

    script:
    """
    ${params.gatk} --java-options '-Xmx${task.memory.giga}G -Xms4G' \
        VariantsToTable \
        -V ${filtered_ann_vcf} \
        -F CHROM \
        -F POS \
        -F REF \
        -F ALT \
        -F HOM-VAR \
        -F HET \
        -F HOM-REF \
        -F NO-CALL \
        -F TYPE \
        -F EVENTLENGTH \
        -F MULTI-ALLELIC \
        -F ANN \
        -GF GT \
        -GF AD \
        -O ${filtered_ann_vcf.baseName}.txt \
        --tmp-dir ${params.scratch_directory}
    """

    stub:
    """
    touch ${filtered_ann_vcf.baseName}.txt
    touch filtered_snps_ann.txt
    touch filtered_indels_ann.txt
    """
}

process VcfToTable {
    // Our own way to turn a VCF file into a R-readable table

    label 'io_intensive'
    debug params.debug

    // Publishing the final resulting .txt files
    publishDir (
        path: "${params.results_directory}/",
        mode: 'copy',
        pattern: "*.txt"
    )

    // Publishing logs from the command execution
    publishDir (
        path: "${params.results_directory}/logs/${task.process}/",
        mode: 'copy',
        pattern: ".command.*"
    )

    // Input: VCF file
    input:
    path(filtered_vcf)

    // Output: TXT tables
    output:
    path("${filtered_vcf.baseName}.txt"), emit: txt
    path(".command.*"), emit: logs

    script:
    """
    vcf_to_table.py \
        --vcf ${filtered_vcf} \
        --output ${filtered_vcf.baseName}.txt \
        --num_allow_missing 0
    """

    stub:
    """
    touch ${filtered_vcf.baseName}.txt
    """
}

// Defining workflow ===========================================================
workflow {

    // Defining samples channel (contains all samples regardless of completion)
    samples = Channel.fromPath(params.metadata)
    .splitCsv(header:true)
    .map { row ->
        fastq1_path = params.samples_directory + row.fastq1
        fastq2_path = row.fastq2 ? params.samples_directory + row.fastq2 : null
        vcf_path = params.vcf_directory +
            row.population + "/" +
            row.population + "_haplotype_caller_results.g.vcf.gz"

        meta = row.subMap(
            'flow_cell',
            'lane',
            'population',
            'barcode',
            'sequencing_facility',
            'internal_library_name'
            )

        def files = [file(fastq1_path, checkIfExists: false)]
        if (fastq2_path) { files << file(fastq2_path, checkIfExists: false) }

        [row.population, meta, files, vcf_path]
    }
    .branch { // branching out into two "subchannels"
        present_vcf: file(it[3]).exists()   // vcf for this file exists
        missing_vcf: !file(it[3]).exists()  // vcf for this file does not exist
    }

    // MAPPING
    // Run mapping processes with missing_vcf samples
    BwaMem( samples.missing_vcf )
    MergeSamFiles(BwaMem.out.bam.groupTuple())
    MarkDuplicates(MergeSamFiles.out.bam)
    BaseRecalibrator(MarkDuplicates.out.bam)
    ApplyBQSR(BaseRecalibrator.out.bam)
    HaplotypeCaller(ApplyBQSR.out.bam)

    // GENOTYPING
    // Extract vcf and tbi from samples.present_vcf
    samples.present_vcf
    .map { sample ->
        def (population, meta, fastq_files, vcf_file) = sample
        def tbi_file = vcf_file + ".tbi"
        return [vcf_file, tbi_file]
    }
    .unique()
    .set { present_vcf_tbi }

    // Mix both channels together. Now we have a full list of our vcfs
    vcf = present_vcf_tbi
    .map { it[0] }
    .mix(HaplotypeCaller.out.vcf)
    .unique()

    tbi = present_vcf_tbi
    .map { it[1] }
    .mix(HaplotypeCaller.out.tbi)
    .unique()

    // Variant Discovery processes
    CombineGVCFs(vcf.collect(), tbi.collect())
    GenotypeGVCFs(CombineGVCFs.out.vcf)
    SelectVariants(GenotypeGVCFs.out.vcf)
    VariantFiltration(SelectVariants.out.vcf)
   
    // // Annotation
    // SnpEff(VariantFiltration.out.vcf.flatten())
    // VariantsToTable(SnpEff.out.vcf)
   
    // Create the filtered_snps.txt file
    VcfToTable(VariantFiltration.out.vcf.flatten())
}