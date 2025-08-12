#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cnv.py
Calculate CNV from coverage data.
This script processes coverage data, calculates log2 ratios against a reference sample,
and visualizes the results using a strip plot.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb


def cnv_calculator(
    data_path: str, reference_sample: str, chromosome: str = None # type: ignore
) -> pd.DataFrame:
    """
    Calculate CNV from coverage data.

    This function reads in corrected coverage data from a TSV file, calculates log2 ratios
    against a reference sample, and visualizes the results using a strip plot.

    Parameters
    ----------
    data_path : str
        Path to the TSV file containing the corrected coverage data.
    reference_sample : str
        Sample ID of the reference sample to use for normalization.

    Returns
    -------
    pd.DataFrame
        DataFrame of the log2(sample coverage / reference coverage).
    """
    # Read in the data
    coverage_data = pd.read_csv(data_path, sep="\t")

    # Pivot the data
    pivoted_data = pd.pivot_table(
        coverage_data,
        index="Sampleid",
        columns="Chromosome",
        values="Corrected Coverage",
    )

    # Calculate the log2 ratio for each sample
    log2_ratios = pivoted_data.copy()
    for chrom in pivoted_data.columns:
        log2_ratios[chrom] = np.log2(
            pivoted_data[chrom] / pivoted_data.at[reference_sample, chrom]
        )

    # If chromosome is given, filter to just that column
    if chromosome:
        if chromosome not in log2_ratios.columns:
            raise ValueError(f"Chromosome '{chromosome}' not found in data.")
        log2_ratios = log2_ratios[[chromosome]]

    # Plot the data using a strip plot
    plt.figure(figsize=(10, 6))
    ax = sb.stripplot(data=log2_ratios, jitter=True, size=5)
    # Add labels for values >= 1
    for i, sample_id in enumerate(log2_ratios.index):
        for j, chromosome in enumerate(log2_ratios.columns):
            value = log2_ratios.iloc[i, j]
            if value >= 1:  # type: ignore
                ax.text(
                    j,  # x position (chromosome index)
                    value + 0.05,  # y position just above the point # type: ignore
                    sample_id,
                    fontsize=8,
                    ha="center"
                )

    plt.axhline(y=0, color="red", linestyle="--")
    plt.axhline(y=1, color="red", linestyle=":")
    plt.xticks(rotation=90)
    plt.ylabel("Log2(Read count ratios)")
    plt.title("Log2(Read count ratios)")

    return log2_ratios


def main():
    """Calculate CNV from coverage data.

    Parameters
    ----------
    data_path : str
        Path to the TSV file containing the corrected coverage data.
    reference_sample : str
        Sample ID of the reference sample to use for normalization.
    output_path : str
        Path to save the CNV results.
    plot_results : bool
        Whether to plot the CNV results.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame containing the log2 ratio of the coverage values for
        each sample relative to the reference sample.
    """
    parser = argparse.ArgumentParser(description="Calculate CNV from coverage data.")
    parser.add_argument(
        "--data-path",
        type=str,
        default="amplicon_covearge.tsv",
        help="Path to the TSV file containing the corrected coverage data.",
    )
    parser.add_argument(
        "--reference-sample",
        type=str,
        default="reference_sample",
        help="Sample ID of the reference sample to use for normalization.",
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default="cnv_results.tsv",
        help="Path to save the CNV results.",
    )
    parser.add_argument(
        "--plot-results",
        action="store_true",
        help="Whether to plot the CNV results.",
    )
    parser.add_argument(
        "--chromosome", 
        type=str,
        help="Specific chromosome to plot")
    args = parser.parse_args()

    # Load the data
    corrected_coverage = pd.read_csv(args.data_path, sep="\t")

    # Ensure the reference sample is in the data
    if args.reference_sample not in corrected_coverage["Sampleid"].values:
        raise ValueError(
            f"Reference sample '{args.reference_sample}' not found in data."
        )
    # Check if the data contains any NaN values
    if corrected_coverage.isna().any().any():
        corrected_coverage.fillna(0, inplace=True)

    # Calculate CNV and plot the results
    cnv_results = cnv_calculator(
        args.data_path,
        args.reference_sample,
        args.chromosome,
    )

    cnv_results.to_csv(args.output_path, sep="\t", index=True, header=True)

    if args.plot_results:
        plt.savefig("cnv_results.png", dpi=300, bbox_inches="tight")

if __name__ == "__main__":
    main()
