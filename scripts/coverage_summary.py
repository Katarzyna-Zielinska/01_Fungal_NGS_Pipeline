
#!/usr/bin/env python3
"""Validate and summarize BEDTools gene coverage for six C. auris samples."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / "results" / "summary"

COVERAGE_FILES = {
    "SRR10059687": ROOT / "results/bedtools/SRR10059687_gene_coverage.tsv",
    "SRR10059688": ROOT / "results/bedtools/SRR10059688_gene_coverage.tsv",
    "SRR10059689": ROOT / "results/bedtools/SRR10059689_gene_coverage.tsv",
    "SRR10059690": ROOT / "results/resistant_cohort/coverage/SRR10059690_gene_coverage.tsv",
    "SRR10059691": ROOT / "results/resistant_cohort/coverage/SRR10059691_gene_coverage.tsv",
    "SRR10059692": ROOT / "results/resistant_cohort/coverage/SRR10059692_gene_coverage.tsv",
}

GROUPS = {
    "SRR10059687": "susceptible",
    "SRR10059688": "susceptible",
    "SRR10059689": "susceptible",
    "SRR10059690": "resistant",
    "SRR10059691": "resistant",
    "SRR10059692": "resistant",
}

EXPECTED_COLUMNS = 8
EXPECTED_GENES = 5583


def load_coverage_file(path: Path, sample: str) -> pd.DataFrame:
    """Load a coverage file and validate its structure."""

    if not path.is_file():
        raise FileNotFoundError(f"Missing coverage file: {path}")

    df = pd.read_csv(path, sep="\t", header=None)

    if df.shape[1] != EXPECTED_COLUMNS:
        raise ValueError(
            f"{sample}: expected {EXPECTED_COLUMNS} columns, "
            f"found {df.shape[1]}"
        )

    if len(df) != EXPECTED_GENES:
        raise ValueError(
            f"{sample}: expected {EXPECTED_GENES} rows, found {len(df)}"
        )

    start = pd.to_numeric(df.iloc[:, 1], errors="raise")
    end = pd.to_numeric(df.iloc[:, 2], errors="raise")
    read_count = pd.to_numeric(df.iloc[:, 4], errors="raise")
    covered_bases = pd.to_numeric(df.iloc[:, 5], errors="raise")
    gene_length = pd.to_numeric(df.iloc[:, 6], errors="raise")
    fraction_covered = pd.to_numeric(df.iloc[:, 7], errors="raise")

    if (end <= start).any():
        raise ValueError(f"{sample}: invalid genomic coordinates")

    # The observed files use column 7 (index 6) as gene length.
    expected_length = end - start
    if not np.array_equal(gene_length.to_numpy(), expected_length.to_numpy()):
        raise ValueError(
            f"{sample}: column 7 does not match BED coordinate lengths"
        )

    if (read_count < 0).any():
        raise ValueError(f"{sample}: negative read count detected")

    if ((covered_bases < 0) | (covered_bases > gene_length)).any():
        raise ValueError(f"{sample}: invalid covered-base count")

    if ((fraction_covered < 0) | (fraction_covered > 1)).any():
        raise ValueError(f"{sample}: coverage fraction outside [0, 1]")

    calculated_fraction = covered_bases / gene_length
    if not np.allclose(
        fraction_covered.to_numpy(),
        calculated_fraction.to_numpy(),
        atol=1e-6,
    ):
        raise ValueError(
            f"{sample}: fraction covered does not match covered bases / length"
        )

    gene_ids = df.iloc[:, 3].astype(str).str.extract(
        r"locus_tag=([^;]+)"
    )[0]

    if gene_ids.isna().any():
        raise ValueError(f"{sample}: missing locus_tag identifiers")

    if gene_ids.duplicated().any():
        raise ValueError(f"{sample}: duplicate gene identifiers")

    df["gene_id"] = gene_ids
    df["read_count"] = read_count
    df["covered_bases"] = covered_bases
    df["gene_length"] = gene_length
    df["fraction_covered"] = fraction_covered

    return df


def summarize_coverage(df: pd.DataFrame, sample: str) -> dict:
    """Calculate descriptive statistics for one sample."""

    read_count = df["read_count"]
    covered_bases = df["covered_bases"]
    gene_length = df["gene_length"]
    fraction_covered = df["fraction_covered"]

    return {
        "sample": sample,
        "group": GROUPS[sample],
        "gene_count": len(df),
        "mean_overlapping_reads_per_gene": float(np.mean(read_count)),
        "median_overlapping_reads_per_gene": float(np.median(read_count)),
        "mean_gene_length_bp": float(np.mean(gene_length)),
        "mean_covered_bases_per_gene": float(np.mean(covered_bases)),
        "median_covered_bases_per_gene": float(np.median(covered_bases)),
        "zero_covered_genes": int((covered_bases == 0).sum()),
        "mean_fraction_covered": float(np.mean(fraction_covered)),
        "genes_with_full_coverage": int((fraction_covered >= 1.0).sum()),
    }


def main() -> None:
    """Validate all samples and save the summary."""

    records = []
    reference_gene_ids = None

    for sample, path in COVERAGE_FILES.items():
        df = load_coverage_file(path, sample)
        current_gene_ids = set(df["gene_id"])

        if reference_gene_ids is None:
            reference_gene_ids = current_gene_ids
        elif current_gene_ids != reference_gene_ids:
            missing = reference_gene_ids - current_gene_ids
            extra = current_gene_ids - reference_gene_ids
            raise ValueError(
                f"{sample}: gene set mismatch; "
                f"missing={len(missing)}, extra={len(extra)}"
            )

        records.append(summarize_coverage(df, sample))
        print(
            f"Validated {sample}: {len(df)} genes; "
            f"{len(current_gene_ids)} unique IDs"
        )

    summary = pd.DataFrame(records).sort_values(
        ["group", "sample"]
    ).reset_index(drop=True)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "coverage_summary.csv"
    summary.to_csv(output_path, index=False)

    print("\nGene coverage summary")
    print(summary.to_string(index=False))

    print("\nGroup summary")
    group_summary = (
        summary.groupby("group")["mean_fraction_covered"]
        .agg(["mean", "std"])
        .reset_index()
    )
    print(group_summary.to_string(index=False))

    print(f"\nAll {len(summary)} samples passed validation.")
    print(f"Shared genes: {len(reference_gene_ids)}")
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()

