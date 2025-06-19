import pandas as pd
import argparse

def convert_vcf_to_csv(vcf_file, output_csv):
    #Read the VCF file, skipping metadata lines starting with '##'
    with open(vcf_file, "r") as file:
        lines = [line.strip() for line in file if not line.startswith("##")]

    #Extract the header and data
    header = lines[0].split("\t")
    data = [line.split("\t") for line in lines[1:]]

    #Create a DataFrame from VCF data
    vcf_df = pd.DataFrame(data, columns=header)

    #Extract the required columns
    position = vcf_df["POS"].astype(int)
    mutants = vcf_df.iloc[:, 9:]  # Columns starting from the 9th are mutants
    mutants.columns = [col + "_sorted" for col in mutants.columns]  # Add "_sorted" to mutant names

    #Create the desired CSV format
    allele_matrix = pd.concat([position, mutants], axis=1)
    allele_matrix.rename(columns={"POS": "Position"}, inplace=True)

    #Save the CSV file
    allele_matrix.to_csv(output_csv, index=False)
    print(f"Converted VCF to allele matrix CSV: {output_csv}")

if __name__ == "__main__":
    # Command-line argument parser
    parser = argparse.ArgumentParser(description="Convert VCF to allele matrix CSV.")
    parser.add_argument("vcf_file", help="Input VCF file")
    parser.add_argument("output_csv", help="Output CSV file")
    args = parser.parse_args()

    # Run the conversion function
    convert_vcf_to_csv(args.vcf_file, args.output_csv)
