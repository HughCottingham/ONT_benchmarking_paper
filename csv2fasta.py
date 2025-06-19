import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse

def csv_to_fasta(input_csv, output_fasta, verbose=False):
    """Convert a CSV allele matrix to a multi-FASTA file."""
    if verbose:
        print(f"Reading input file: {input_csv}")

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(input_csv)

    # Initialize a dictionary to store sequences for each sample
    sequences = {col: [] for col in df.columns[1:]}  # Skip the "Position" column

    # Concatenate alleles for each sample
    for index, row in df.iterrows():
        for sample in sequences.keys():
            sequences[sample].append(row[sample])

    # Convert sequences to FASTA records
    records = [
        SeqRecord(Seq("".join(alleles)), id=sample, description="")
        for sample, alleles in sequences.items()
    ]

    # Write the records to a FASTA file
    SeqIO.write(records, output_fasta, "fasta")

    if verbose:
        print(f"FASTA file written to {output_fasta}")

# Set up command-line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert a CSV allele matrix to a multi-FASTA file.")
    parser.add_argument("input_csv", help="Path to the input CSV file.")
    parser.add_argument("output_fasta", help="Path to the output FASTA file.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Print verbose output.")

    args = parser.parse_args()

    # Run the conversion
    csv_to_fasta(args.input_csv, args.output_fasta, args.verbose)

