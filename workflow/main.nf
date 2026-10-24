
nextflow.enable.dsl=2

params.raw_dir = 'data/raw/SRR26185932'
params.outdir = 'results/nextflow'
params.index_prefix = 'data/reference/hisat2_index_current/candida_auris_B8441_v2'
params.genes_bed = 'results/bedtools/candida_auris_genes.sorted.bed'


process FASTQC_RAW {

    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("fastqc_raw")

    script:
    """
    mkdir -p fastqc_raw
    fastqc $read1 $read2 --outdir fastqc_raw
    """

    stub:
    """
    mkdir -p fastqc_raw
    touch fastqc_raw/${sample}_1_fastqc.html
    touch fastqc_raw/${sample}_2_fastqc.html
    """
}


process CUTADAPT {

    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2), path(qc_dir)

    output:
    tuple val(sample),
          path("${sample}_1.cutadapt.fastq"),
          path("${sample}_2.cutadapt.fastq")

    script:
    """
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o ${sample}_1.cutadapt.fastq \
        -p ${sample}_2.cutadapt.fastq \
        $read1 $read2 \
        --minimum-length 20 \
        --cores 2
    """

    stub:
    """
    touch ${sample}_1.cutadapt.fastq
    touch ${sample}_2.cutadapt.fastq
    """
}


process FASTP {

    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_1.fastq"),
          path("${sample}_2.fastq")

    script:
    """
    fastp \
        -i $read1 \
        -I $read2 \
        -o ${sample}_1.fastq \
        -O ${sample}_2.fastq \
        --disable_adapter_trimming \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 20 \
        --n_base_limit 5 \
        --length_required 50 \
        --thread 2 \
        --html ${sample}_fastp.html \
        --json ${sample}_fastp.json
    """

    stub:
    """
    touch ${sample}_1.fastq
    touch ${sample}_2.fastq
    touch ${sample}_fastp.html
    touch ${sample}_fastp.json
    """
}


process FASTQC_TRIMMED {

    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample), path("fastqc_trimmed")

    script:
    """
    mkdir -p fastqc_trimmed
    fastqc $read1 $read2 --outdir fastqc_trimmed
    """

    stub:
    """
    mkdir -p fastqc_trimmed
    touch fastqc_trimmed/${sample}_1_fastqc.html
    touch fastqc_trimmed/${sample}_2_fastqc.html
    """
}


process HISAT2 {

    tag "$sample"

    input:
    tuple val(sample), path(read1), path(read2), path(qc_dir)

    output:
    tuple val(sample),
          path("${sample}.sam"),
          path("${sample}.hisat2.summary.txt")

    script:
    """
    hisat2 \
        -x ${params.index_prefix} \
        -1 $read1 \
        -2 $read2 \
        -S ${sample}.sam \
        --summary-file ${sample}.hisat2.summary.txt \
        -p 4
    """

    stub:
    """
    touch ${sample}.sam
    echo "Stub run: HISAT2 was not executed" \
        > ${sample}.hisat2.summary.txt
    """
}


process SAMTOOLS {

    tag "$sample"

    input:
    tuple val(sample), path(sam), path(hisat_summary)

    output:
    tuple val(sample),
          path("${sample}.sorted.bam"),
          path("${sample}.sorted.bam.bai"),
          path("${sample}.flagstat.txt"),
          path("${sample}.stats.txt")

    script:
    """
    samtools sort \
        -@ 2 \
        -o ${sample}.sorted.bam \
        $sam

    samtools index \
        ${sample}.sorted.bam

    samtools flagstat \
        ${sample}.sorted.bam \
        > ${sample}.flagstat.txt

    samtools stats \
        ${sample}.sorted.bam \
        > ${sample}.stats.txt
    """

    stub:
    """
    touch ${sample}.sorted.bam
    touch ${sample}.sorted.bam.bai
    touch ${sample}.flagstat.txt
    touch ${sample}.stats.txt
    """
}


process BEDTOOLS_COVERAGE {

    tag "$sample"

    cpus 1
    memory '2 GB'

    publishDir "${params.outdir}/bedtools", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(bam), path(bai)
    path genes_bed

    output:
    path "${sample}.gene_coverage.tsv"

    script:
    """
    set -euo pipefail

    test -s "$bam"
    test -s "$bai"
    test -s "$genes_bed"
    samtools quickcheck "$bam"

    OUT="${sample}.gene_coverage.tsv"
    : > "\$OUT"

    samtools idxstats "$bam" \
        | awk '\$1 != "*" {print \$1}' \
        | while read -r ref; do

            echo "Przetwarzam referencję: \$ref" >&2

            tmp_bed="\${ref}.genes.bed"
            tmp_bam="\${ref}.bam"

            awk -F '\\t' -v r="\$ref" '\$1 == r' "$genes_bed" > "\$tmp_bed"

            if [ ! -s "\$tmp_bed" ]; then
                rm -f "\$tmp_bed"
                continue
            fi

            samtools view -b -o "\$tmp_bam" "$bam" "\$ref"
            samtools quickcheck "\$tmp_bam"

            bedtools coverage \
                -a "\$tmp_bed" \
                -b "\$tmp_bam" \
                >> "\$OUT"

            rm -f "\$tmp_bed" "\$tmp_bam"
        done

    test -s "\$OUT"
    test "\$(wc -l < "\$OUT")" -eq "\$(wc -l < "$genes_bed")"
    """

    stub:
    """
    printf "stub\\n" > "${sample}.gene_coverage.tsv"
    """
}


workflow {

    log.info "Project 1 — Fungal NGS Pipeline"
    log.info "Raw data: ${params.raw_dir}"
    log.info "Output: ${params.outdir}"
    log.info "HISAT2 index: ${params.index_prefix}"
    log.info "Gene annotation: ${params.genes_bed}"

    raw_reads = Channel.of(
        tuple(
            'SRR26185932',
            file("${params.raw_dir}/SRR26185932_1.fastq"),
            file("${params.raw_dir}/SRR26185932_2.fastq")
        )
    )

    genes_bed = file(params.genes_bed)

    FASTQC_RAW(raw_reads)

    raw_with_qc = raw_reads.join(FASTQC_RAW.out)

    cutadapt_reads = CUTADAPT(raw_with_qc)

    fastp_reads = FASTP(cutadapt_reads)

    FASTQC_TRIMMED(fastp_reads)

    trimmed_with_qc = fastp_reads.join(FASTQC_TRIMMED.out)

    HISAT2(trimmed_with_qc)

    SAMTOOLS(HISAT2.out)

    BEDTOOLS_COVERAGE(SAMTOOLS.out.map { sample, bam, bai, flagstat, stats ->
        tuple(sample, bam, bai)
    }, genes_bed)
}
