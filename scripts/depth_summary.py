
#!/usr/bin/env python3
"""Estimate genome-wide sequencing depth from SAMtools stats reports."""

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
STATS_DIR = ROOT / "results" / "bam_qc"
OUTPUT_DIR = ROOT / "results" / "summary"

GENOME_FILE = ROOT / "data" / "reference" / "candida_auris.genome"

SAMPLES = {
    "SRR10059687": "susceptible",
    "SRR10059688": "susceptible",
    "SRR10059689": "susceptible",
    "SRR10059690": "resistant",
    "SRR10059691": "resistant",
    "SRR10059692": "resistant",
}


def calculate_genome_length(path: Path) -> tuple[int, int]:
    """Return the number of contigs and total reference length in bp."""

    if not path.is_file():
        raise FileNotFoundError(f"Missing genome length file: {path}")

    contigs = 0
    genome_length = 0

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()

            if not fields:
                continue

            if len(fields) < 2:
                raise ValueError(f"Invalid genome file line: {line.rstrip()}")

            length = int(fields[1])

            if length <= 0:
                raise ValueError(f"Invalid contig length: {line.rstrip()}")

            contigs += 1
            genome_length += length

    if genome_length == 0:
        raise ValueError("Genome length is zero.")

    return contigs, genome_length


def extract_stat(text: str, label: str) -> int:
    """Extract an integer value from a SAMtools SN statistics line."""

    pattern = rf"^SN\s+{re.escape(label)}:\s+(\d+)"
    match = re.search(pattern, text, flags=re.MULTILINE)

    if match is None:
        raise ValueError(f"Statistic not found: {label}")

    return int(match.group(1))


def main() -> None:
    contigs, genome_length = calculate_genome_length(GENOME_FILE)
    records = []

    for sample, group in SAMPLES.items():
        stats_path = STATS_DIR / f"{sample}_stats.txt"

        if not stats_path.is_file():
            raise FileNotFoundError(
                f"Missing SAMtools stats report: {stats_path}"
            )

        text = stats_path.read_text(encoding="utf-8")

        mapped_bases = extract_stat(text, "bases mapped (cigar)")
        mapped_reads = extract_stat(text, "reads mapped")
        primary_reads = extract_stat(text, "raw total sequences")

        estimated_depth = mapped_bases / genome_length

        records.append({
            "sample": sample,
            "group": group,
            "reference_contigs": contigs,
            "genome_length_bp": genome_length,
            "primary_reads": primary_reads,
            "mapped_reads": mapped_reads,
            "mapped_bases_cigar": mapped_bases,
            "estimated_mean_depth_x": estimated_depth,
        })

    df = pd.DataFrame(records)
    df = df.sort_values(["group", "sample"]).reset_index(drop=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "depth_summary.csv"
    df.to_csv(output_path, index=False)

    print("Estimated genome-wide sequencing depth")
    print(f"Reference contigs: {contigs}")
    print(f"Reference length: {genome_length:,} bp")
    print()
    print(df.to_string(index=False, formatters={
        "estimated_mean_depth_x": lambda value: f"{value:.2f}x"
    }))

    print("\nSummary by group")
    group_summary = (
        df.groupby("group")["estimated_mean_depth_x"]
        .agg(["mean", "std", "min", "max"])
        .reset_index()
    )
    print(group_summary.to_string(index=False, formatters={
        column: lambda value: f"{value:.2f}x"
        for column in ["mean", "std", "min", "max"]
    }))

    print("\nNote: this is an estimated mean depth, not per-base coverage.")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
