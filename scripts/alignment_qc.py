#!/usr/bin/env python3
"""Summarize SAMtools flagstat metrics for the Candida auris cohort."""

from pathlib import Path
import re

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
QC_DIR = PROJECT_ROOT / "results" / "bam_qc"
OUTPUT_DIR = PROJECT_ROOT / "results" / "summary"

SAMPLES = {
    "SRR10059687": "susceptible",
    "SRR10059688": "susceptible",
    "SRR10059689": "susceptible",
    "SRR10059690": "resistant",
    "SRR10059691": "resistant",
    "SRR10059692": "resistant",
}


def extract_metric(text: str, pattern: str) -> tuple[int, float | None]:
    """Extract a read count and optional percentage from flagstat text."""
    match = re.search(pattern, text, flags=re.MULTILINE)
    if not match:
        raise ValueError(f"Metric not found: {pattern}")

    count = int(match.group(1))
    percentage = float(match.group(2)) if match.lastindex and match.lastindex >= 2 else None
    return count, percentage


def parse_flagstat(path: Path, sample: str, group: str) -> dict:
    """Parse selected metrics from one samtools flagstat report."""
    text = path.read_text(encoding="utf-8")

    total, _ = extract_metric(
        text, r"^(\d+) \+ \d+ in total"
    )
    primary, _ = extract_metric(
        text, r"^(\d+) \+ \d+ primary$"
    )
    primary_mapped, mapping_rate = extract_metric(
        text, r"^(\d+) \+ \d+ primary mapped \(([\d.]+)%"
    )
    properly_paired, paired_rate = extract_metric(
        text, r"^(\d+) \+ \d+ properly paired \(([\d.]+)%"
    )

    return {
        "sample": sample,
        "group": group,
        "total_records": total,
        "primary_records": primary,
        "primary_mapped": primary_mapped,
        "primary_mapping_rate_pct": mapping_rate,
        "properly_paired": properly_paired,
        "properly_paired_pct": paired_rate,
    }


def main() -> None:
    """Build the cohort-level QC table."""
    records = []

    for sample, group in SAMPLES.items():
        flagstat_path = QC_DIR / f"{sample}_flagstat.txt"
        if not flagstat_path.is_file():
            raise FileNotFoundError(f"Missing flagstat report: {flagstat_path}")

        records.append(parse_flagstat(flagstat_path, sample, group))

    df = pd.DataFrame(records)
    df = df.sort_values(["group", "sample"]).reset_index(drop=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_csv = OUTPUT_DIR / "alignment_qc_summary.csv"
    df.to_csv(output_csv, index=False)

    print("Alignment QC summary")
    print(df.to_string(index=False))
    print()
    print(f"Mean primary mapping rate: "
          f"{np.mean(df['primary_mapping_rate_pct']):.2f}%")
    print(f"SD primary mapping rate: "
          f"{np.std(df['primary_mapping_rate_pct'], ddof=1):.2f}%")
    print(f"Saved: {output_csv}")


if __name__ == "__main__":
    main()
