nextflow.enable.dsl=2

params.outdir = 'results/resistant_cohort'
params.index_prefix = 'data/reference/hisat2_index_current/candida_auris_B8441_v2'
params.genes_bed = 'results/bedtools/candida_auris_genes.bam_order.bed'
params.genome_file = 'data/reference/candida_auris.genome'

process SRA_TO_FASTQ {
    tag "$sample"
    maxForks 1
    cpus 4

    publishDir "${params.outdir}/fastq/raw", mode: "copy"

    input:
    tuple val(sample), path(sra)

    output:
    tuple val(sample), path("${sample}_1.fastq"), path("${sample}_2.fastq"), emit: reads


    script:
    """
    set -euo pipefail
    mkdir -p tmp_sra

    fasterq-dump "$sra" \
        --split-files \
        --threads ${task.cpus} \
        --outdir . \
        --temp tmp_sra

    test -s "${sample}_1.fastq"
    test -s "${sample}_2.fastq"
    """

    stub:
    """
    touch "${sample}_1.fastq" "${sample}_2.fastq"
    """
}

process FASTQC_RAW {
    tag "$sample"
    cpus 2

    publishDir "${params.outdir}/fastqc/raw", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("${sample}_*_fastqc.html"), emit: html
    path("${sample}_*_fastqc.zip"), emit: zip


    script:
    """
    fastqc "$read1" "$read2" --threads ${task.cpus}
    """

    stub:
    """
    touch "${sample}_1_fastqc.html" "${sample}_2_fastqc.html"
    touch "${sample}_1_fastqc.zip" "${sample}_2_fastqc.zip"
    """
}

process CUTADAPT {
    tag "$sample"
    cpus 2

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("${sample}_1.cutadapt.fastq"), path("${sample}_2.cutadapt.fastq"), emit: reads

    script:
    """
    set -euo pipefail

    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "${sample}_1.cutadapt.fastq" \
        -p "${sample}_2.cutadapt.fastq" \
        "$read1" "$read2" \
        --minimum-length 20 \
        --cores ${task.cpus} \
        > "${sample}.cutadapt.log" 2>&1

    test -s "${sample}_1.cutadapt.fastq"
    test -s "${sample}_2.cutadapt.fastq"
    """

    stub:
    """
    touch "${sample}_1.cutadapt.fastq" "${sample}_2.cutadapt.fastq"
    """
}

process FASTP {
    tag "$sample"
    cpus 2

    publishDir "${params.outdir}/fastq/trimmed", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_1.fastq"),
          path("${sample}_2.fastq"),
          path("${sample}_fastp.html"),
          path("${sample}_fastp.json"),
          emit: reads


    script:
    """
    set -euo pipefail

    fastp \
        -i "$read1" \
        -I "$read2" \
        -o "${sample}_1.fastq" \
        -O "${sample}_2.fastq" \
        --disable_adapter_trimming \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 20 \
        --n_base_limit 5 \
        --length_required 50 \
        --thread ${task.cpus} \
        --html "${sample}_fastp.html" \
        --json "${sample}_fastp.json"

    test -s "${sample}_1.fastq"
    test -s "${sample}_2.fastq"
    """

    stub:
    """
    touch "${sample}_1.fastq" "${sample}_2.fastq"
    touch "${sample}_fastp.html" "${sample}_fastp.json"
    """
}

process FASTQC_TRIMMED {
    tag "$sample"
    cpus 2

    publishDir "${params.outdir}/fastqc/trimmed", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path("${sample}_*_fastqc.html"), emit: html
    path("${sample}_*_fastqc.zip"), emit: zip


    script:
    """
    fastqc "$read1" "$read2" --threads ${task.cpus}
    """

    stub:
    """
    touch "${sample}_1_fastqc.html" "${sample}_2_fastqc.html"
    touch "${sample}_1_fastqc.zip" "${sample}_2_fastqc.zip"
    """
}

process HISAT2_SORT {
    tag "$sample"
    maxForks 1
    cpus 4
    memory '8 GB'

    publishDir "${params.outdir}/alignment", mode: "copy"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai"),
          path("${sample}.hisat2.summary.txt"),
          path("${sample}.flagstat.txt"),
          emit: alignment


    script:
    """
    set -euo pipefail

    hisat2 \
        -x "${params.index_prefix}" \
        -1 "$read1" \
        -2 "$read2" \
        --summary-file "${sample}.hisat2.summary.txt" \
        -p ${task.cpus} \
    | samtools sort -@ 2 -o "${sample}.sorted.bam" -

    samtools index "${sample}.sorted.bam"
    samtools quickcheck -v "${sample}.sorted.bam"
    samtools flagstat "${sample}.sorted.bam" > "${sample}.flagstat.txt"

    test -s "${sample}.sorted.bam"
    test -s "${sample}.sorted.bam.bai"
    """

    stub:
    """
    touch "${sample}.sorted.bam" "${sample}.sorted.bam.bai"
    touch "${sample}.hisat2.summary.txt" "${sample}.flagstat.txt"
    """
}

process BEDTOOLS_COVERAGE {
    tag "$sample"
    maxForks 1
    cpus 1
    memory '4 GB'

    publishDir "${params.outdir}/bedtools", mode: "copy"

    input:
    tuple val(sample), path(bam), path(bai), path(summary), path(flagstat)
    path genes_bed
    path genome_file

    output:
    tuple val(sample), path("${sample}_gene_coverage.tsv"), emit: coverage


    script:
    """
    set -euo pipefail

    samtools quickcheck -v "$bam"
    test -s "$bai"
    test -s "$genes_bed"

    OUT="${sample}_gene_coverage.tsv"
    : > "\$OUT"

    bedtools coverage \
        -sorted \
        -g "$genome_file" \
        -a "$genes_bed" \
        -b "$bam" \
        > "\$OUT"

    cut -f1-4 "$genes_bed" | sort -k1,1 -k2,2n -k3,3n -k4,4 > expected.coords
    cut -f1-4 "\$OUT" | sort -k1,1 -k2,2n -k3,3n -k4,4 > observed.coords
    diff -q expected.coords observed.coords
    """

    stub:
    """
    touch "${sample}_gene_coverage.tsv"
    """
}

workflow {
    raw_sra = Channel.of(
        tuple('SRR10059690', file('data/raw/SRR10059690/SRR10059690.sra')),
        tuple('SRR10059691', file('data/raw/SRR10059691/SRR10059691.sra')),
        tuple('SRR10059692', file('data/raw/SRR10059692/SRR10059692.sra'))
    )

    genes_bed = file(params.genes_bed)
    genome_file = file(params.genome_file)

    SRA_TO_FASTQ(raw_sra)

    FASTQC_RAW(SRA_TO_FASTQ.out.reads)

    CUTADAPT(SRA_TO_FASTQ.out.reads)

    FASTP(CUTADAPT.out.reads)

    FASTQC_TRIMMED(FASTP.out.reads.map {
        sample, r1, r2, html, json -> tuple(sample, r1, r2)
    })

    HISAT2_SORT(FASTP.out.reads.map {
        sample, r1, r2, html, json -> tuple(sample, r1, r2)
    })

    BEDTOOLS_COVERAGE(HISAT2_SORT.out.alignment, genes_bed, genome_file)
}
