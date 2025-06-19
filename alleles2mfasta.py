import argparse
from Bio import SeqIO
import csv

def create_mfasta(reference_fasta, allele_matrix, output_mfasta):
    #Load the reference sequence
    reference_seq = None
    with open(reference_fasta, "r") as ref_file:
        for record in SeqIO.parse(ref_file, "fasta"):
            reference_seq = list(record.seq)
            reference_id = record.id
            break

    #Parse the SNP allele matrix
    mutants = []
    positions = []

    with open(allele_matrix, "r") as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # Read header row
        mutants = header[2:]  # Mutant names
        positions = []

        # Prepare SNP data for each mutant
        snp_data = {mutant: {} for mutant in mutants}

        for row in reader:
            pos = int(row[0]) - 1  # Convert 1-based position to 0-based index
            ref_base = row[1]

            for mutant, allele in zip(mutants, row[2:]):
                if allele != ref_base:  # Only record if allele differs from the reference
                    snp_data[mutant][pos] = allele

    #Create modified sequences for each mutant
    mutant_sequences = {}
    for mutant in mutants:
        mutant_seq = reference_seq[:]  # Copy reference sequence
        for pos, allele in snp_data[mutant].items():
            mutant_seq[pos] = allele  # Apply SNP
        mutant_sequences[mutant] = "".join(mutant_seq)

    #Write the multi-FASTA file
    with open(output_mfasta, "w") as out_fasta:
        for mutant, seq in mutant_sequences.items():
            out_fasta.write(f">{mutant}\n")
            out_fasta.write(f"{seq}\n")

    print(f"Multi-FASTA file created: {output_mfasta}")

if __name__ == "__main__":
    # Set up command-line arguments
    parser = argparse.ArgumentParser(description="Create a multi-FASTA file with SNP profiles.")
    parser.add_argument("-r", "--reference", required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("-a", "--alleles", required=True, help="Path to the SNP allele matrix (CSV format).")
    parser.add_argument("-o", "--output", required=True, help="Output multi-FASTA file.")

    args = parser.parse_args()

    # Run the function
    create_mfasta(args.reference, args.alleles, args.output)
