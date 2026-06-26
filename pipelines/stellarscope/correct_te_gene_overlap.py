#!/usr/bin/env python3

"""
Correct overlap between STARsolo GeneFull assignments and Stellarscope TE counts.

Purpose:
    Identify molecules that were counted both as genes (STARsolo)
    and as transposable elements (Stellarscope) using Cell Barcode (CB)
    and UMI (UB), and generate a correction table to remove potential
    double counting.

Input:
    - STARsolo BAM file with CB/UB tags
    - Stellarscope TE BAM file with CB/UB tags
    - STARsolo GeneFull count matrix
    - Stellarscope TE count matrix

Output:
    - Correction table containing TE counts to subtract
"""


import pysam
import pandas as pd
import argparse


def extract_star_molecules(star_bam):
    """
    Extract molecules counted by STARsolo using CB + UB identifiers.
    """

    molecules = set()

    bam = pysam.AlignmentFile(star_bam, "rb")

    for read in bam:

        if read.has_tag("CB") and read.has_tag("UB"):

            cb = read.get_tag("CB")
            ub = read.get_tag("UB")

            molecules.add(f"{cb}_{ub}")

    bam.close()

    return molecules



def find_te_overlap(te_bam, star_molecules):
    """
    Identify Stellarscope molecules that were already assigned
    to genes by STARsolo.
    """

    overlap = []

    bam = pysam.AlignmentFile(te_bam, "rb")

    for read in bam:

        if read.has_tag("CB") and read.has_tag("UB"):

            cb = read.get_tag("CB")
            ub = read.get_tag("UB")

            molecule_id = f"{cb}_{ub}"


            if molecule_id in star_molecules:

                overlap.append({
                    "Cell": cb,
                    "UMI": ub
                })


    bam.close()

    return pd.DataFrame(overlap)



def main():

    parser = argparse.ArgumentParser(
        description="Remove STARsolo/Stellarscope molecule overlap"
    )


    parser.add_argument(
        "--star_bam",
        required=True,
        help="STARsolo BAM file"
    )

    parser.add_argument(
        "--te_bam",
        required=True,
        help="Stellarscope TE BAM file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output correction table"
    )


    args = parser.parse_args()


    print("Extracting STARsolo molecules...")

    star_molecules = extract_star_molecules(
        args.star_bam
    )


    print(
        f"STARsolo molecules detected: {len(star_molecules)}"
    )


    print("Searching TE overlap...")


    overlap = find_te_overlap(
        args.te_bam,
        star_molecules
    )


    overlap.to_csv(
        args.output,
        sep="\t",
        index=False
    )


    print(
        f"Saved correction table: {args.output}"
    )



if __name__ == "__main__":
    main()
