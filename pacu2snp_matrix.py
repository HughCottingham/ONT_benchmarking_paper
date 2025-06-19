import argparse
import pandas as pd

def generate_allele_matrix(snp_matrix_file, snp_positions_file, output_file):
    #Parse the SNP positions file
    positions = pd.read_csv(snp_positions_file, sep="\t")
    positions_list = positions['position'].tolist()

    #Parse the SNP matrix FASTA file
    mutant_data = {}
    with open(snp_matrix_file, "r") as fasta:
        current_mutant = None
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                current_mutant = line[1:]  # Extract mutant name
                mutant_data[current_mutant] = []
            else:
                mutant_data[current_mutant].extend(list(line))  # Store sequence

    #Create the allele matrix
    allele_matrix = pd.DataFrame({'Position': positions_list})
    for mutant, bases in mutant_data.items():
        allele_matrix[mutant] = bases

    #Save the matrix to a CSV file
    allele_matrix.to_csv(output_file, index=False)
    print(f"Allele matrix saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an allele matrix from SNP matrix and positions.")
    parser.add_argument("-m", "--matrix", required=True, help="Path to the SNP matrix FASTA file.")
    parser.add_argument("-p", "--positions", required=True, help="Path to the SNP positions TSV file.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for the allele matrix.")
    args = parser.parse_args()

    generate_allele_matrix(args.matrix, args.positions, args.output)
