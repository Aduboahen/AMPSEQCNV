#!/usr/bin/env python
"""
gc_bias_correction.py
Correct GC bias in read coverage following Benjamini & Speed (2012), 
and summarise per-amplicon statistics.

Usage:
    python gc_bias_correction.py --bam sample.bam --fasta reference.fa --bed targets.bed \
          --window 300 --output corrected_coverage.tsv --summary amplicon_summary.tsv

Dependencies:
    pysam
    pyfaidx
    numpy
    statsmodels
    scipy
    bedtools (for preparing BED file externally if needed)
"""

import argparse
import csv
import sys
from collections import defaultdict
import pysam
from pyfaidx import Fasta
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d


def compute_gc(seq: str) -> float:
    """
    Compute GC content of a sequence as a percentage.

    Args:
        seq (str): A DNA sequence string.

    Returns:
        float: GC content as a percentage. Returns NaN if the sequence is empty.
    """
    # seq = seq
    if len(seq) == 0:
        return np.nan
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def load_bed_regions(bed_path: str):
    """
    Load BED regions from a file and return as a list of (chrom, start, end) tuples.

    Args:
        bed_path (str): Path to the BED file.

    Returns:
        list: A list of (chrom, start, end) tuples.
    """
    regions = []
    with open(bed_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or line.strip() == "":
                continue
            chrom, start, end = line.strip().split()[:3]
            regions.append((chrom, int(start), int(end)))
    return regions


def collect_coverage_in_bed(
    bam_path: str,
    fa: Fasta,
    regions,
    window: int,
):
    """
    Collect GC content and coverage within BED regions from a BAM file.

    Args:
        bam_path (str): Path to the input BAM file.
        fa (Fasta): A Fasta object for the reference genome.
        regions (list): A list of (chrom, start, end) tuples for the regions of interest.
        window (int): The window size for computing GC content.
        min_mapq (int): The minimum mapping quality to include. Defaults to 0.
        max_depth (int): The maximum read depth to include. Defaults to 100000.

    Returns:
        tuple: A tuple of four arrays: GC content, coverage, positions, and region indices.
    """

    bam = pysam.AlignmentFile(bam_path, "rb")  # type: ignore
    gc_vals = []
    cov_vals = []
    positions = []
    region_map = []

    sys.stderr.write("Collecting coverage and GC content within BED regions...\n")
    for region_index, (chrom, start, end) in enumerate(regions):
        for pileup in bam.pileup(chrom, start, end, stepper="all"):
            pos = pileup.pos  # type: ignore
            cov = pileup.nsegments
            half = window // 2
            seq_start = max(0, pos - half)
            seq_end = pos + half  # + window % 2
            try:
                seq = (fa.get_seq(chrom, seq_start, seq_end)).seq  # type: ignore
            except KeyError:
                continue
            gc = compute_gc(seq)
            if np.isnan(gc):
                continue
            gc_vals.append(gc)
            cov_vals.append(cov)
            positions.append((chrom, pos))
            region_map.append(region_index)
    bam.close()
    return np.array(gc_vals), np.array(cov_vals), positions, region_map


def fit_loess(gc_vals: np.ndarray, cov_vals: np.ndarray, frac: float = 0.3):
    """
    Fit a LOESS curve to the given GC content and coverage arrays.

    Args:
        gc_vals (ndarray): GC content values.
        cov_vals (ndarray): Coverage values.
        frac (float): Fraction of data to use for LOESS regression. Defaults to 0.3.

    Returns:
        ndarray: The LOESS curve as a 2D array of (GC, coverage) values.
    """
    sys.stderr.write("Fitting LOESS curve...\n")
    loess = sm.nonparametric.lowess(cov_vals, gc_vals, frac=frac, return_sorted=True)
    return loess


def correct_coverage(
    gc_vals: np.ndarray,
    cov_vals: np.ndarray,
    loess_result,
    output_path: str,
    positions,
    region_map,
    region_names,
    summary_path=None,
):
    """
    Write the corrected coverage and expected coverage to a TSV file.

    Parameters
    ----------
    gc_vals : ndarray
        GC content values.
    cov_vals : ndarray
        Coverage values.
    loess_result : ndarray
        LOESS curve as a 2D array of (GC, coverage) values.
    output_path : str
        Path to write the corrected coverage output file.
    positions : list
        List of (chrom, pos) tuples.
    region_map : list
        List of region indices for each position.
    region_names : list
        List of region names.
    summary_path : str, optional
        Path to write the per-amplicon summary output file. Defaults to None.
    """
    x, y = loess_result[:, 0], loess_result[:, 1]
    interp = interp1d(x, y, bounds_error=False, fill_value="extrapolate")  # type: ignore
    sys.stderr.write(f"Writing corrected coverage to {output_path} ...\n")

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["Chrom", "Pos", "GC", "Coverage", "Expected", "Corrected", "Region"]
        )
        for gc, cov, (chrom, pos), region_idx in zip(
            gc_vals, cov_vals, positions, region_map
        ):
            expected = interp(gc)
            corrected = (
                cov / expected if expected != 0 and not np.isnan(expected) else 0
            )
            region_name = (
                region_names[region_idx]
                if region_idx < len(region_names)
                else f"region_{region_idx}"
            )
            writer.writerow(
                [
                    chrom,
                    pos,
                    f"{gc:.2f}",
                    cov,
                    f"{expected:.4f}",
                    f"{corrected:.4f}",
                    region_name,
                ]
            )

    if summary_path:
        sys.stderr.write(f"Writing per-amplicon summary to {summary_path} ...\n")
        per_region_cov = defaultdict(list)
        per_region_exp = defaultdict(list)
        for corr, exp, region_idx in zip(cov_vals, interp(gc_vals), region_map):
            per_region_cov[region_idx].append(corr)
            per_region_exp[region_idx].append(exp)

        with open(summary_path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(
                [
                    "Chromosome",
                    "Region",
                    "Coverage",
                    "Expected Coverage",
                    "Corrected Coverage",
                    # "NumBases",
                ]
            )
            for (idx, values), (values2) in zip(
                per_region_cov.items(), per_region_exp.values()
            ):
                region_name = (
                    region_names[idx] if idx < len(region_names) else f"region_{idx}"
                )
                writer.writerow(
                    [
                        region_name.split(":")[0],
                        region_name.split(":")[1],
                        f"{np.nansum(values):.4f}",
                        f"{np.nansum(values2):.4f}",
                        f"{np.nansum(values)/np.nansum(values2, ):.4f}",
                        # len(values),
                    ]
                )


def main():
    """
    GC bias correction and per-amplicon coverage summary.

    This script takes a BAM file of aligned reads, a FASTA file of the reference genome,
    and a BED file of target amplicon regions. It then computes GC content and read coverage
    within the regions, fits a LOESS curve to the data, and applies the LOESS curve to
    correct the coverage values. The corrected coverage values are then written to a TSV
    file. Optionally, a per-amplicon summary of mean and median corrected coverage values
    can be written to a separate TSV file.

    """
    parser = argparse.ArgumentParser(
        description="GC bias correction and per-amplicon coverage summary."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file (sorted).")
    parser.add_argument(
        "--fasta",
        required=True,
        help="Reference genome in FASTA format (indexed by faidx).",
    )
    parser.add_argument(
        "--bed", required=True, help="BED file of target amplicon regions."
    )
    parser.add_argument(
        "--window",
        type=int,
        default=300,
        help="Fragment window size to compute GC (%(default)s).",
    )
    parser.add_argument(
        "--frac",
        type=float,
        default=0.3,
        help="LOESS smoothing parameter (%(default)s).",
    )
    parser.add_argument(
        "--output", default="corrected_coverage.tsv", help="Output TSV path."
    )
    parser.add_argument(
        "--summary",
        default="per_amplicon.tsv",
        help="Output TSV path for per-amplicon summary.",
    )
    parser.add_argument(
        "--sample",
        help="Sample name for output files.",
    )
    # parser.add_argument(
    # "--min-mapq", type=int, default=0, help="Minimum mapping quality to include."
    # )
    args = parser.parse_args()

    fa = Fasta(args.fasta, sequence_always_upper=True)
    regions = load_bed_regions(args.bed)
    region_names = [f"{c}:{s}-{e}" for c, s, e in regions]

    gc_vals, cov_vals, positions, region_map = collect_coverage_in_bed(
        args.bam, fa, regions, args.window
    )

    if len(gc_vals) == 0:
        sys.stderr.write("No positions collected; check inputs.\n")
        return

    loess_result = fit_loess(gc_vals, cov_vals, args.frac)
    correct_coverage(
        gc_vals,
        cov_vals,
        loess_result,
        (
            f"{args.sample}_corrected_coverage_per_site.tsv"
            if args.sample
            else args.output
        ),
        positions,
        region_map,
        region_names,
        (
            f"{args.sample}_corrected_coverage_per_amplicon.tsv"
            if args.sample
            else args.summary
        ),
    )

    sys.stderr.write("Done.\n")


if __name__ == "__main__":
    main()
