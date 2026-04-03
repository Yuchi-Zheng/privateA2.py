import gzip
import os
import argparse
from collections import defaultdict, Counter
import numpy as np

def parse_population_file(pop_file):
    """Parse population assignment file, return a dictionary mapping individuals to populations"""
    indv_to_pop = {}
    with open(pop_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            indv = parts[0]
            pop = parts[1]
            indv_to_pop[indv] = pop
    return indv_to_pop

def count_private_alleles(vcf_path, indv_to_pop):
    """Count the number of private alleles for each population and for each individual"""
    # Initialize data structures
    pop_allele_counts = defaultdict(lambda: defaultdict(int))
    private_allele_counts = defaultdict(int)
    pop_samples = defaultdict(list)
    # New: Individual private allele counts
    individual_private_counts = defaultdict(int)
    # New: For storing private allele information per site
    site_private_info = []
    
    # Create mapping from population to samples
    for indv, pop in indv_to_pop.items():
        pop_samples[pop].append(indv)
    
    # Open VCF file
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    with open_func(vcf_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith('##'):
                continue  # Skip metadata
            
            # Process header line
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                all_samples = header[9:]
                
                # Ensure all samples in VCF are defined in the population file
                for sample in all_samples:
                    if sample not in indv_to_pop:
                        print(f"Warning: Sample '{sample}' not defined in population file, will be skipped")
                
                # Create mapping from sample index to population
                sample_to_pop = {}
                sample_idx_to_name = {}  # New: Mapping from index to sample name
                for idx, sample in enumerate(all_samples):
                    if sample in indv_to_pop:
                        sample_to_pop[idx] = indv_to_pop[sample]
                        sample_idx_to_name[idx] = sample  # New: Store mapping from index to sample name
                continue
            
            # Process variant site line
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alts = fields[4].split(',')
            format_fields = fields[8].split(':')
            gt_index = format_fields.index('GT') if 'GT' in format_fields else 0
            alleles = [ref] + alts
            
            # Initialize allele occurrence for this site
            allele_in_pop = {i: set() for i in range(len(alleles))}
            # New: Record genotype for each sample at this site
            sample_genotypes = {}
            
            # Process each sample's genotype
            for idx, sample_data in enumerate(fields[9:]):
                if idx not in sample_to_pop:
                    continue  # Skip samples without population definition
                
                pop = sample_to_pop[idx]
                sample_name = sample_idx_to_name[idx]  # Get sample name
                sample_fields = sample_data.split(':')
                gt_str = sample_fields[gt_index]
                
                # Parse genotype (handle various formats)
                if '/' in gt_str:
                    alleles_idx = gt_str.split('/')
                elif '|' in gt_str:
                    alleles_idx = gt_str.split('|')
                else:
                    alleles_idx = [gt_str]
                
                # Process alleles in genotype
                called_alleles = []
                for a in alleles_idx:
                    if a == '.' or a == './.' or a == '.|.':
                        continue
                    try:
                        allele_idx = int(a)
                        if 0 <= allele_idx < len(alleles):
                            called_alleles.append(allele_idx)
                            allele_in_pop[allele_idx].add(pop)
                    except ValueError:
                        continue
                
                # Record genotype for this sample
                sample_genotypes[sample_name] = called_alleles
            
            # Check each allele if it's a private allele
            private_alleles_in_site = []
            for allele_idx, pops in allele_in_pop.items():
                if len(pops) == 1:  # Only appears in one population
                    pop_name = next(iter(pops))
                    private_allele_counts[pop_name] += 1
                    private_alleles_in_site.append((allele_idx, pop_name))
            
            # New: If this site has private alleles, update individual counts
            if private_alleles_in_site:
                # For each private allele, find its population
                for allele_idx, pop_name in private_alleles_in_site:
                    # Check which individuals in this population carry this allele
                    for sample_name, called_alleles in sample_genotypes.items():
                        sample_pop = indv_to_pop.get(sample_name)
                        if sample_pop == pop_name and allele_idx in called_alleles:
                            # This individual carries this private allele, count +1
                            individual_private_counts[sample_name] += 1
    
    return private_allele_counts, pop_samples, individual_private_counts  # Modified: Return three values

def main():
    parser = argparse.ArgumentParser(description='Calculate the number of private alleles for each population and for each individual')
    parser.add_argument('--vcf', required=True, help='Input VCF file path (supports .gz compressed format)')
    parser.add_argument('--pop', required=True, help='Population assignment file path')
    parser.add_argument('--output', default='private_alleles.csv', help='Output CSV file path')
    parser.add_argument('--indiv-output', default='individual_private_alleles.csv', 
                       help='Individual private alleles output file path (default: individual_private_alleles.csv)')
    args = parser.parse_args()
    
    print(f"Parsing population file: {args.pop}")
    indv_to_pop = parse_population_file(args.pop)
    print(f"Found {len(indv_to_pop)} individuals assigned to {len(set(indv_to_pop.values()))} populations")
    
    print(f"Processing VCF file: {args.vcf}")
    private_alleles, pop_samples, individual_private_counts = count_private_alleles(args.vcf, indv_to_pop)  # Modified: Receive three return values
    
    # Calculate sample count for each population
    pop_sample_counts = {pop: len(samples) for pop, samples in pop_samples.items()}
    
    # Write population-level CSV results
    with open(args.output, 'w') as f:
        f.write("Population,SampleCount,PrivateAlleles\n")
        for pop in sorted(private_alleles.keys()):
            sample_count = pop_sample_counts.get(pop, 0)
            alleles_count = private_alleles[pop]
            f.write(f"{pop},{sample_count},{alleles_count}\n")
    
    # New: Write individual-level CSV results
    with open(args.indiv_output, 'w') as f:
        f.write("Individual,Population,PrivateAlleles\n")
        # Sort by population and individual name
        sorted_individuals = sorted(individual_private_counts.keys(), 
                                   key=lambda x: (indv_to_pop.get(x, ""), x))
        for indv in sorted_individuals:
            pop = indv_to_pop.get(indv, "Unknown")
            alleles_count = individual_private_counts[indv]
            f.write(f"{indv},{pop},{alleles_count}\n")
    
    print(f"Population-level results saved to: {args.output}")
    print(f"Individual-level results saved to: {args.indiv_output}")
    
    print("\nPopulation Private Allele Statistics:")
    print(f"{'Population':<15}{'SampleCount':<10}{'PrivateAlleles':<15}")
    for pop in sorted(private_alleles.keys()):
        print(f"{pop:<15}{pop_sample_counts.get(pop,0):<10}{private_alleles[pop]:<15}")
    
    # New: Output individual statistics summary
    print("\nIndividual Private Allele Statistics Summary:")
    if individual_private_counts:
        # Group statistics by population
        pop_indiv_counts = defaultdict(list)
        for indv, count in individual_private_counts.items():
            pop = indv_to_pop.get(indv, "Unknown")
            pop_indiv_counts[pop].append(count)
        
        print(f"{'Population':<15}{'Individuals':<10}{'AvgPrivateAlleles':<20}{'Range':<15}")
        for pop in sorted(pop_indiv_counts.keys()):
            counts = pop_indiv_counts[pop]
            if counts:
                avg_count = sum(counts) / len(counts)
                min_count = min(counts)
                max_count = max(counts)
                print(f"{pop:<15}{len(counts):<10}{avg_count:<20.2f}{f'{min_count}-{max_count}':<15}")
    else:
        print("No individual private alleles found")

if __name__ == "__main__":
    main()