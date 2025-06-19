import pandas as pd
import re
import argparse
import io
import subprocess
import os

def load_depth_data(depth_file):
    if os.path.exists(depth_file):
        return pd.read_csv(depth_file, delimiter='\t', header=None, names=['#CHROM', 'POS', 'DEPTH'])
    else:
        print("Warning: Depth file {} does not exist.".format(depth_file))
        return pd.DataFrame(columns=['#CHROM', 'POS', 'DEPTH'])

def replace_values(df, regex_pattern, depth_threshold):
    pattern = str(regex_pattern)
    depth_data = {}
    for col in df.columns:
        if col not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
            depth_file = "{}_depth_filtered.txt".format(col)
            depth_data[col] = load_depth_data(depth_file)
    
    rows_to_drop = []
    total_columns = len(df.columns) - 9  # Total relevant columns
    total_rows = len(df.index)
    quartiles = [total_rows // 4, total_rows // 2, 3 * total_rows // 4, total_rows]
    rows_iterated = 0

    for row_idx in df.index:
        low_depth_count = 0
        for col in df.columns:
            if col not in ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']:
                if col in depth_data:
                    depth_df = depth_data[col]
                    depth_row = depth_df[(depth_df['#CHROM'] == df.at[row_idx, '#CHROM']) &
                                         (depth_df['POS'] == df.at[row_idx, 'POS'])]
                    
                    if not depth_row.empty:
                        depth = depth_row.iloc[0]['DEPTH']
                        if depth < depth_threshold:
                            low_depth_count += 1
                            df.at[row_idx, col] = 'N'
                            continue
                cell_value = str(df.at[row_idx, col])
                if str(pattern) in cell_value:
                    df.at[row_idx, col] = df.at[row_idx, 'REF']
                else:
                    alt_values = df.at[row_idx, 'ALT'].split(',')
                    first_integer = int(re.search(r'\d', cell_value).group())
                    df.at[row_idx, col] = alt_values[first_integer - 1]

        if low_depth_count > 0.05 * total_columns:
            rows_to_drop.append(row_idx)
        
        rows_iterated += 1
        if int(rows_iterated) in quartiles:
            print("Processed {} variant calls out of {}. Variant calls filtered due to <95 percent conservation across isolates: {}".format(rows_iterated, total_rows, len(rows_to_drop)))

    df.drop(rows_to_drop, inplace=True)
    return df

def filter_rows_by_pattern(df):
    """
    Removes rows where any column contains a pattern matching ':<float between 0.4 and 0.6>'
    with up to 4 decimal places, followed by a tab separator.
    Prints a message for each filtered row.
    """
    pattern = r":0\.([4-5]\d{0,3}|6{0,4})\t"
    rows_to_keep = []

    for idx, row in df.iterrows():
        row_string = row.to_csv(sep='\t', index=False)  # Convert row to tab-separated string
        match = re.search(pattern, row_string)
        if match:
            print(f"Filtering out row {idx} due to match: {match.group(0)}")
        else:
            rows_to_keep.append(idx)

    return df.loc[rows_to_keep]

def write_fasta(df, fasta_file):
    with open(fasta_file, 'w') as file:
        for col in df.columns:
            if col not in ['POS', 'REF']:
                file.write(">{}\n".format(col))
                file.write("".join(df[col].astype(str)) + "\n")

def run_snp_sites(input_fasta, output_fasta):
    result = subprocess.run(['snp-sites', '-c', input_fasta], stdout=open(output_fasta, 'w'), stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode != 0:
        print("Error running snp-sites: {}".format(result.stderr))
        raise subprocess.CalledProcessError(result.returncode, 'snp-sites')

def run_snp_dists(fasta_file, distance_matrix_file):
    with open(distance_matrix_file, 'w') as output_file:
        result = subprocess.run(['snp-dists', fasta_file], stdout=output_file, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            print("Error running snp-dists: {}".format(result.stderr))
            raise subprocess.CalledProcessError(result.returncode, 'snp-dists')

def run_fasttree(fasta_file, tree_file):
    with open(tree_file, 'w') as output_file:
        result = subprocess.run(['Fasttree', '-gtr', fasta_file], stdout=output_file, stderr=subprocess.PIPE, universal_newlines=True)
        if result.returncode != 0:
            print("Error running Fasttree: {}".format(result.stderr))
            raise subprocess.CalledProcessError(result.returncode, 'Fasttree')

def run_gubbins(fasta_file, prefix):
    print("Running gubbins...")
    result = subprocess.run(['run_gubbins.py', '-p', prefix, '-f 50', fasta_file], stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode != 0:
        print("Error running Gubbins: {}".format(result.stderr))
        raise subprocess.CalledProcessError(result.returncode, 'run_gubbins.py')

def filter_lines(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    filtered_lines = [line for line in lines if not line.startswith('##')]
    return filtered_lines

parser = argparse.ArgumentParser(description='Process a table and replace values based on regex pattern.')
parser.add_argument('input_file', type=str, help='Path to the input TSV file')
parser.add_argument('regex_pattern', type=str, help='Regex pattern to match')
parser.add_argument('output_file', type=str, help='Path to the output TSV file')
parser.add_argument('--chrom_match', type=str, help='String to match in the #CHROM column', default=None)
parser.add_argument('--depth_threshold', type=int, help='Depth value below which to issue a warning', default=10)

args = parser.parse_args()

filtered_lines = filter_lines(args.input_file)
df = pd.read_csv(io.StringIO(''.join(filtered_lines)), delimiter='\t')

df = df[~df['FILTER'].str.contains('RefCall', na=False)]

if args.chrom_match:
    df = df[df['#CHROM'] == args.chrom_match]

df = replace_values(df, args.regex_pattern, args.depth_threshold)
df = filter_rows_by_pattern(df)
df = df.drop(columns=['#CHROM', 'ID', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'])
df.to_csv(args.output_file, sep=',', index=False)
