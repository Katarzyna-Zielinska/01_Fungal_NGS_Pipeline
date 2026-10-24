
nextflow.enable.dsl=2

params.sample = 'SRR26185932'
params.bam = 'data/bam/SRR26185932.sorted.bam'
params.bai = 'data/bam/SRR26185932.sorted.bam.bai'
params.bed = 'results/bedtools/candida_auris_genes.sorted.bed'
params.outdir = 'results/nextflow/bedtools_verified'

process BEDTOOLS_COVERAGE {

    tag "${sample}"

    cpus 1
    memory '2 GB'

    publishDir params.outdir, mode: 'copy', overwrite: false

    input:
    val sample
    path bam
    path bai
    path bed

    output:
    path "${sample}.gene_coverage.tsv", emit: coverage

    script:
    """
    set -euo pipefail

    test -s "$bam"
    test -s "$bai"
    test -s "$bed"
    samtools quickcheck "$bam"

    OUT="${sample}.gene_coverage.tsv"
    : > "\$OUT"

    samtools idxstats "$bam" \
        | awk '\$1 != "*" {print \$1}' \
        | while read -r ref; do

            echo "Przetwarzam referencję: \$ref" >&2

            tmp_bed="\${ref}.genes.bed"
            tmp_bam="\${ref}.bam"

            awk -F '\\t' -v r="\$ref" '\$1 == r' "$bed" > "\$tmp_bed"

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
    test "\$(wc -l < "\$OUT")" -eq "\$(wc -l < "$bed")"
    """

    stub:
    """
    printf "stub\\n" > "${sample}.gene_coverage.tsv"
    """
}

workflow {
    BEDTOOLS_COVERAGE(
        params.sample,
        file(params.bam),
        file(params.bai),
        file(params.bed)
    )
}
