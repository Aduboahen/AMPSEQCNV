#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gc_bias_correction.py
Correct GC bias in read coverage following Benjamini & Speed (2012), 
and summarise per-amplicon statistics.

Usage:
    python gc_bias_correction.py --bam BAMS --fasta reference.fa --bed targets.bed \
          --window 300 --output corrected_coverage.tsv --summary amplicon_summary.tsv

Dependencies:
    pysam
    pyfaidx
    numpy
    statsmodels
    scipy
"""
import os
import argparse
import csv
import sys
from collections import defaultdict
from pysam import AlignmentFile
from pyfaidx import Fasta
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d


def compute_gc(seq: str) -> float:
    """
    Calculate the GC content of a DNA sequence.

    This function computes the percentage of guanine (G) and cytosine (C)
    bases in a given DNA sequence.

    Parameters
    ----------
    seq : str
        The DNA sequence for which to calculate GC content.

    Returns
    -------
    float
        The GC content as a percentage of the total sequence length.
        Returns NaN if the sequence is empty.
    """

    if not seq:
        return np.nan
    gc = seq.count("G") + seq.count("C")
    return 100.0 * gc / len(seq)


def load_bed_regions(bed_path: str):
    """
    Load regions from a BED file.

    Parameters
    ----------
    bed_path : str
        Path to the BED file.

    Returns
    -------
    list
        A list of 3-tuples containing the chromosome, start position, and
        end position of each region in the BED file.
    """
    regions = []
    try:
        with open(bed_path, encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("#") or line.strip() == "":
                    continue
                chrom, start, end = line.strip().split()[:3]
                regions.append((chrom, int(start), int(end)))
    except FileNotFoundError:
        sys.stderr.write(f"Error: BED file {bed_path} not found.\n")
        sys.exit(1)
    return regions


def collect_coverage_in_bed(bam_path: str, fa: Fasta, regions, window: int):
    """
    Collect coverage and GC content within regions from a BAM file.

    Parameters
    ----------
    bam_path : str
        Path to the BAM file.
    fa : Fasta
        Reference genome sequence.
    regions : list
        List of 3-tuples containing the chromosome, start position, and
        end position of each region.
    window : int
        Window size of the fragment to compute GC content.

    Returns
    -------
    tuple
        A tuple of three numpy arrays: the GC content, coverage values, and
        region mapping of each position.
    """
    bam = AlignmentFile(bam_path, "rb")
    gc_vals = []
    cov_vals = []
    positions = []
    region_map = []

    sys.stderr.write("Collecting coverage and GC content within BED regions...\n")
    for region_index, (chrom, start, end) in enumerate(regions):
        for pileup in bam.pileup(chrom, start, end, stepper="all"):
            pos = pileup.pos  # pyright: ignore[reportAttributeAccessIssue]
            cov = pileup.nsegments
            half = window // 2
            seq_start = max(0, pos - half)
            seq_end = pos + half
            try:
                seq = fa.get_seq(
                    chrom, seq_start, seq_end
                ).seq  # pyright: ignore[reportAttributeAccessIssue]
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
    return np.array(gc_vals), np.array(cov_vals), region_map


def fit_loess(gc_vals: np.ndarray, cov_vals: np.ndarray, frac: float = 0.3):
    """
    Fit a LOESS curve to the coverage vs GC content data.

    Parameters
    ----------
    gc_vals : numpy.ndarray
        A 1D array of GC content values.
    cov_vals : numpy.ndarray
        A 1D array of coverage values.
    frac : float, optional
        The fraction of the data to use when computing the LOESS
        curve. Default is 0.3.

    Returns
    -------
    tuple
        A tuple of two 1D numpy arrays: the x and y values of the
        LOESS curve.
    """
    sys.stderr.write("Fitting LOESS curve...\n")
    loess = sm.nonparametric.lowess(cov_vals, gc_vals, frac=frac, return_sorted=True)

    return loess


def correct_coverage(
    gc_vals: np.ndarray,
    cov_vals: np.ndarray,
    loess_result,
    region_map,
    region_names,
    sampleid: str,
    summary_path: str = "amplicon_coverage.tsv",
):
    """
    Correct coverage values using a LOESS curve.

    Parameters
    ----------
    gc_vals : numpy.ndarray
        A 1D array of GC content values.
    cov_vals : numpy.ndarray
        A 1D array of coverage values.
    loess_result : tuple
        A tuple of two 1D numpy arrays: the x and y values of the
        LOESS curve.
    region_map : list
        A list of region indices that the coverage values correspond
        to.
    region_names : list
        A list of region names (chromosome:start-end) that the
        coverage values correspond to.
    sampleid : str
        The sample ID to write to the output file.
    summary_path : str, optional
        The path to the output file to write the corrected coverage
        values to. If not provided, does not write to a file.

    """
    x, y = loess_result[:, 0], loess_result[:, 1]
    interp = interp1d(
        x,
        y,
        bounds_error=False,
        fill_value="extrapolate",  # pyright: ignore[reportArgumentType]
    )

    per_region_cov = defaultdict(list)
    per_region_exp = defaultdict(list)

    # Clip predicted expected coverage to >= 0
    expected_covs = np.clip(interp(gc_vals), 0, None)

    for corr, exp, region_idx in zip(cov_vals, expected_covs, region_map):
        per_region_cov[region_idx].append(corr)
        per_region_exp[region_idx].append(exp)

    if summary_path:
        with open(summary_path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            for region_idx in sorted(per_region_cov.keys()):
                values = per_region_cov[region_idx]
                values2 = per_region_exp[region_idx]
                denom = np.nansum(values2)
                ratio = 0 if denom == 0 else np.nansum(values) / denom
                region_name = (
                    region_names[region_idx]
                    if region_idx < len(region_names)
                    else f"region_{region_idx}"
                )
                writer.writerow(
                    [
                        region_name.split(":")[0],
                        region_name.split(":")[1],
                        f"{np.nansum(values):.2f}",
                        f"{denom:.2f}",
                        f"{ratio:.2f}" if not np.isnan(ratio) else 0,
                        sampleid,
                    ]
                )


def get_output_filename(path: str) -> str:
    """Check if file exists, prompt to delete or rename."""
    while os.path.exists(path):
        choice = (
            input(f"File '{path}' already exists. [D]elete, [R]ename, or [Q]uit? ")
            .strip()
            .lower()
        )

        if choice == "d":
            os.remove(path)
            print(f"Deleted '{path}'.")
        elif choice == "r":
            path = input("Enter new filename: ").strip()
        elif choice == "q":
            raise SystemExit("Aborted by user.")
        else:
            print("Invalid choice. Please enter D, R, or Q.")

    return path


def main():
    """
    Entry point for GC bias correction and per-amplicon coverage summary.

    Collects coverage and GC content within regions from a BAM file, 
    fits a LOESS curve to the data, and corrects the coverage values. 
    Computes per-amplicon coverage summary and writes to an output file.

    Parameters
    ----------
    bams : str
        Path to sorted BAM files.
    fasta : str
        Reference genome in FASTA format (indexed by faidx).
    bed : str
        BED file of target amplicon regions.
    window : int, optional
        Fragment window size to compute GC (default is 300).
    frac : float, optional
        LOESS smoothing parameter (default is 0.3).
    summary : str, optional
        Output TSV path for per-amplicon summary (default is "amplicon_coverage.tsv").

    Returns
    -------
    None
    """
    parser = argparse.ArgumentParser(
        description="GC bias correction and per-amplicon coverage summary."
    )
    parser.add_argument(
        "--bams", required=True, help="Path to directory containing sorted BAM files."
    )
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
        "--summary",
        default="amplicon_coverage.tsv",
        help="Output TSV path for per-amplicon summary.",
    )
    args = parser.parse_args()

    summary_path = get_output_filename(args.summary)

    with open(summary_path, "a", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            [
                "Chromosome",
                "Region",
                "Coverage",
                "Expected Coverage",
                "Corrected Coverage",
                "Sampleid",
            ]
        )
    sys.stderr.write(f"Writing per-amplicon summary to {args.summary} ...\n")

    try:
        fa = Fasta(args.fasta, sequence_always_upper=True)
    except FileNotFoundError:
        sys.stderr.write(f"Error: FASTA file {args.fasta} not found.\n")
        sys.exit(1)

    regions = load_bed_regions(args.bed)
    region_names = [f"{c}:{s}-{e}" for c, s, e in regions]

    for root, _, files in os.walk(args.bams):
        for file in files:
            if file.endswith(".bam"):
                bam = os.path.join(root, file)
                sampleid = file.split(".")[0]
                gc_vals, cov_vals, region_map = collect_coverage_in_bed(
                    bam, fa, regions, args.window
                )
                if len(gc_vals) == 0:
                    sys.stderr.write("No positions collected; check inputs.\n")
                    continue
                loess_result = fit_loess(gc_vals, cov_vals, args.frac)
                correct_coverage(
                    gc_vals,
                    cov_vals,
                    loess_result,
                    region_map,
                    region_names,
                    sampleid,
                    args.summary,
                )

    sys.stderr.write("Done.\n")


if __name__ == "__main__":
    main()
