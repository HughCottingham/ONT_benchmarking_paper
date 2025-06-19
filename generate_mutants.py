import random
from Bio import SeqIO
from Bio.Seq import Seq

# Input files
fasta_file = "reference.fasta"
bed_file = "core_regions.bed"

#Select 1000 random positions in the core genome
def get_positions_from_bed(bed_file):
    positions = []
    with open(bed_file, "r") as bed:
        for line in bed:
            chrom, start, end = line.strip().split()[:3]
            positions.extend(range(int(start), int(end) + 1))
    return random.sample(positions, 1000)

positions = get_positions_from_bed(bed_file)

#Read the FASTA file
def read_fasta(fasta_file):
    record = SeqIO.read(fasta_file, "fasta")
    return record

record = read_fasta(fasta_file)

#Assign reference alleles and generate alternative alleles
def get_reference_bases(record, positions):
    seq = record.seq
    return {pos: seq[pos] for pos in positions}

def generate_alternative_alleles(ref_base):
    possible_bases = [b for b in "ACGT" if b != ref_base]
    return random.choice(possible_bases)

reference_bases = get_reference_bases(record, positions)
alternative_bases = {pos: generate_alternative_alleles(ref) for pos, ref in reference_bases.items()}

#Generate SNP profiles with evenly distributed SNP distances
def generate_snp_profiles(num_mutants, num_positions, max_distance):
    profiles = []
    
    # Create a list of SNP counts that are evenly distributed from 0 to max_distance
    snp_counts = [int(i * max_distance / (num_mutants - 1)) for i in range(num_mutants)]
    
    # Shuffle the SNP counts for randomness
    random.shuffle(snp_counts)
    
    # For each mutant, assign a set of SNPs based on the corresponding SNP count
    snp_sites = random.sample(range(num_positions), num_positions)  # Randomize the order of SNP sites
    for snp_count in snp_counts:
        selected_snp_positions = random.sample(snp_sites, snp_count)
        profiles.append(selected_snp_positions)
    
    return profiles

snp_profiles = generate_snp_profiles(10, len(positions), 600) ## I chose 600 here as the max distance as that was the typical max I saw in real data with ~1000 SNP sites across the whole ST

#Create 10 mutant FASTA files
def mutate_sequence(record, positions, snp_profile, reference_bases, alternative_bases, output_file):
    seq = list(record.seq)  # Convert sequence to a mutable list
    for pos_idx in range(len(positions)):
        pos = positions[pos_idx]
        if pos_idx in snp_profile:
            # Use alternative allele for mutation
            seq[pos] = alternative_bases[pos]
        else:
            # Keep the reference allele
            seq[pos] = reference_bases[pos]
    record.seq = Seq("".join(seq))  # Convert back to a Seq object
    SeqIO.write(record, output_file, "fasta")

for i, snp_profile in enumerate(snp_profiles, 1):
    output_file = f"mutant_{i}.fasta"
    mutate_sequence(record, positions, snp_profile, reference_bases, alternative_bases, output_file)

print("10 mutant FASTA files generated with evenly distributed SNP distances")

