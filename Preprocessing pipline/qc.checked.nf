nextflow.enable.dsl=2

params.ref = ''                         // Reference genome file path
params.reads = ''                       // Pattern to match input FASTQ files
params.indexdir = './index'             // Output directory
params.alignedir = './alignment'        // Output directory
params.countdir = './counts'            // Output directory
params.single_end = false               // Boolean flag for single-end reads
params.annotation = ''                  // Annotation file path

// Process to index the reference genome
process indexReference {
    publishDir "${params.indexdir}", mode: 'copy'

    input:
    path reference

    output:
    path "*.bt2", emit: index_files

    script:
    """
    prefix=\$(basename ${reference} .fasta)
    bowtie2-build ${reference} \$prefix
    """
}

// Process to align paired-end reads using Bowtie2
process alignReads {
    input:
    tuple val(sample_id), path(read1), path(read2)
    path index_files

    output:
    path "${sample_id}.sam"

    script:
    """
    idx_base=\$(basename ${index_files[0]} .1.bt2)
    bowtie2 -x \${idx_base} -1 ${read1} -2 ${read2} -S ${sample_id}.sam
    """
}

// Process to align single-end reads using Bowtie2
process alignReads_single_end {
    input:
    tuple val(sample_id), path(read1)
    path index_files

    output:
    path "${sample_id}.sam"

    script:
    """
    idx_base=\$(basename ${index_files[0]} .1.bt2)
    bowtie2 -x \${idx_base} -U ${read1} -S ${sample_id}.sam
    """
}

// Process to convert SAM to BAM
process samToBam {
    input:
    path sam

    output:
    path "${sam.baseName}.bam"

    script:
    """
    samtools view -Sb ${sam} > ${sam.baseName}.bam
    """
}

params.sort_by_name = params.single_end.toBoolean() ? true : false

process sorting {
    publishDir "${params.alignedir}", mode: 'copy'

    input:
    path bam

    output:
    path "${bam.baseName}_sorted.bam"

    script:
    def sort_flag = params.sort_by_name ? "-n" : ""
    """
    samtools sort ${sort_flag} ${bam} -o ${bam.baseName}_sorted.bam
    """
}

process alignmentStats {
    publishDir "${params.alignedir}/stats", mode: 'copy'

    input:
    path bam

    output:
    path "${bam.simpleName}_flagstat.txt"

    script:
    """
    samtools flagstat ${bam} > ${bam.simpleName}_flagstat.txt
    """
}

process indexBam {
    input:
    path sorted_bam

    output:
    path "${sorted_bam}.bai"

    script:
    """
    samtools index ${sorted_bam}
    """
}

// Process to count features using featureCounts for paired-end reads
process featureCounts {
    publishDir "${params.countdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path annotation

    output:
    path "${sample_id}_counts.txt"
    path "${sample_id}_counts.txt.summary"

    script:
    """
    featureCounts -p --countReadPairs -t gene -g gene_id -a ${annotation} -o ${sample_id}_counts.txt ${bam_file}
    """
}

// Process to count features using featureCounts for single-end reads
process featureCounts_single_end {
    publishDir "${params.countdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)
    path annotation

    output:
    path "${sample_id}_counts.txt"
    path "${sample_id}_counts.txt.summary"

    script:
    """
    featureCounts -t gene -g gene_id -a ${annotation} -o ${sample_id}_counts.txt ${bam_file}
    """
}

// Main workflow
workflow {
    // Validate required parameters
    if (!params.ref) {
        error "Reference genome file path (params.ref) must be provided."
    }
    if (!params.reads) {
        error "Input reads pattern (params.reads) must be provided."
    }
    if (!params.annotation) {
        error "Annotation file (params.annotation) must be provided."
    }

    reference_ch = Channel.fromPath(params.ref)
    annotation_ch = Channel.fromPath(params.annotation)

    ref_index = indexReference(reference_ch)

    if (params.single_end.toBoolean()) {
        reads_ch = Channel
            .fromFilePairs(params.reads, flat: true)
            .map { sample_id, read1 -> tuple(sample_id, read1) }

        sam_ch = alignReads_single_end(reads_ch, ref_index.index_files)
        bam_ch = samToBam(sam_ch)
        sorted_bam_ch = sorting(bam_ch)
        stats_ch = alignmentStats(sorted_bam_ch)

        counts_ch = featureCounts_single_end(
            sorted_bam_ch.map { bam -> tuple(bam.baseName.replace('_sorted.bam',''), bam) },
            annotation_ch
        )
    } else {
        reads_ch = Channel
            .fromFilePairs(params.reads)
            .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }

        sam_ch = alignReads(reads_ch, ref_index.index_files)
        bam_ch = samToBam(sam_ch)
        sorted_bam_ch = sorting(bam_ch)
        stats_ch = alignmentStats(sorted_bam_ch)

        counts_ch = featureCounts(
            sorted_bam_ch.map { bam -> tuple(bam.baseName.replace('_sorted.bam',''), bam) },
            annotation_ch
        )
    }
}

