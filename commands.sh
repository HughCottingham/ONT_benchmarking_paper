### Get read subsets ###

for f in *fastq.gz;do base=$(basename $f .fastq.gz);seqkit fx2tab -n -i -l $f > $base"_read_lengths".txt;done
for f in *lengths.txt;do base=$(basename $f .txt);sort -R $f > $base"_randomised.txt";rm $f;done
for f in *_randomised.txt;do base=$(basename $f _randomised.txt);python3 print_depth_levels.py $f $genome_size $base;done
for fq in *.fastq.gz; do base=$(basename "$fq" .fastq.gz); for txt in ${base}_*x_reads.txt; do depth=$(basename "$txt" | sed -E 's/.*_([0-9]+x)_reads.txt/\1/'); seqtk subseq "$fq" "$txt" | gzip > "${base}_${depth}.fastq.gz"; done; done



### Generating and checking mutant genomes with 1000 SNP sites in the core genome for simulated variant calling ###

#Identify core genome in reference genome based on pangenome annotation across its ST. Current directory includes all hybrid assemblies from a given ST. Reference genome is named reference.fasta.
for assembly in *.fasta; do     basename=$(basename "$assembly" .fasta);     prokka --outdir prokka_output/$basename --cpus 32 --prefix $basename $assembly; done
panaroo -i prokka_output/*/*.gff -o panaroo_output --clean-mode strict --threads 32
num_genomes=13 ## Change to number of genomes in ST
awk -F',' -v genomes=$num_genomes 'NR > 1 {
    count = 0;
    for (i = 4; i <= genomes + 3; i++) if ($i != "") count++;
    if (count / genomes >= 0.95) print $1
}' panaroo_output/gene_presence_absence.csv > core_gene_list.txt
grep -F -f core_gene_list.txt prokka_output/reference/reference.gff > core_regions.gff
awk '$3 == "CDS" {print $1"\t"$4-1"\t"$5}' core_regions.gff > core_regions.bed

#Generate mutant genomes
python generate_mutants.py

#Check their SNP profiles to make sure everything looks good
for mutated in mutant*.fasta; do     minimap2 -ax asm5 reference.fasta $mutated > ${mutated%.fasta}_aligned.sam; done
for f in *sam;do base=$(basename $f _aligned.sam);samtools view -Sb $f | samtools sort -o $base"_sorted.bam"; done
for f in *bam;do base=$(basename $f _sorted.bam);bcftools mpileup -f reference.fasta $f | bcftools call -mv -Oz -o $f".vcf.gz";done
for f in *gz;do echo $f;bcftools view -H $f | grep -v '^#' | wc -l; done




### Variant Calling ###

#Minimap2 alignment to generate sorted bam files - input for PACU and Clair3

for f in *fastq.gz;do base=$(basename $f .fastq.gz);minimap2 -t 16 --cs --MD -aLx map-ont reference.fasta $f > $base".bam";done
for f in *bam;do base=$(basename $f .bam);samtools sort $f > $base"_sorted.bam";done
for f in *sorted.bam;do base=$(basename $f _sorted.bam);rm $base".bam";done
for f in *bam;do samtools index $f;done


#PACU - PACU runs gubbins v3 so no need to run separately like for Clair3 and SKA2. Input directory contains sorted bam files

PACU --ont-in . --ref-fasta reference.fasta --output output --dir-working work --threads 16 --min-snp-qual 10 --min-snp-depth 10 --min-snp-dist 0 --min-snp-af 0.7
python pacu2snpmatrix.py -m output/snp_matrix.fasta -p output/snp_positions.tsv -o pacu_allele_matrix.csv


# Clair3 - use same gubbins version as PACU

for f in *sorted.bam;do base=$(basename $f _sorted.bam);samtools depth -aa $f > $base"_depth.txt";run_clair3.sh --bam_fn=$f --ref_fn=reference.fasta --snp_min_af=0.75 --qual=5 --threads=4 --platform="ont" --model_path=r1041_e82_400bps_sup_v430 --output=$base"_clair3_out" --sample_name=$base --haploid_precise --no_phasing_for_fa;done
for f in *_clair3_out;do bcftools norm $f/merge_output.vcf.gz -f reference.fasta -a -c e -m - |bcftools norm -aD |bcftools +remove-overlaps - |bcftools filter -e 'abs(ILEN)>0 || QUAL<5' -o $base"_bcf_out.vcf";done
for f in *vcf;do bgzip $f
for f in *vcf.gz;do bcftools index $f
bcftools merge *vcf.gz -o all_filtered_snps.vcf
awk '!/^#/ {print $2}' all_filtered_snps.vcf > positions.txt
for f in *depth.txt;do base=$(basename $f _depth.txt);awk 'NR==FNR{positions[$1]; next} $2 in positions' positions.txt $f > $base"_depth_filtered.txt";done
python3 vcf2snpmatrix.py all_filtered_snps.vcf ".:." alleles_cons0.95.csv --chrom_match "contig_1" --depth_threshold 10
python3 alleles2mfasta.py -r reference.fasta -a alleles_cons0.95.csv -o wga.fasta
run_gubbins.py --threads 8 -p wga wga.fasta;done
python3 vcf2csv.py wga.summary_of_snp_distribution.vcf clair3_allele_matrix.csv

# SKA2

paste <(ls *.fasta | sed 's/[.].*$//g') <(ls -d *.fasta) > ska_input.tsv;ska build -f ska_input.tsv -k 31 -o ska_index --threads 4;ska map -o ska_map.aln --ambig-mask --threads 4 reference.fasta ska_index.skf
run_gubbins.py --prefix wga ska_map.aln --threads 32
python3 vcf2csv.py wga.summary_of_snp_distribution.vcf ska2_allele_matrix.csv

# Reddog - ran with default settings then:

python3 alleles2mfasta.py -r reference.fasta -a alleles_cons0.95.csv -o wga.fasta
run_gubbins.py --prefix wga ska_map.aln --threads 32
python3 vcf2csv.py wga.summary_of_snp_distribution.vcf reddog_allele_matrix.csv


# Pairwise SNP distances from matrices

python3 csv2fasta.py tool_allele_matrix.csv tool_msa.fasta
snp-dists tool_msa.fasta > tool_dists.tsv
