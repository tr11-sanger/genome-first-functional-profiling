import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Parse BAM to generate taxonomic and functional profiles. Reassigns multiple mappings with winner-takes-all strategy. Calculates coverage depth and breadth for each species, genome and CDS.')
parser.add_argument('-s', "--sourmash_profile", type=str,
                    required=True,
                    help="Sourmash profile filepath.")
parser.add_argument('-g', "--genome_species", type=str,
                    required=True,
                    help="Genome to species mapping TSV filepath.")
parser.add_argument('-o', "--output_fp", type=str,
                    required=True,
                    help="Output filepath.")
parser.add_argument('-n', "--n_genomes_per_species", type=int,
                    default=4,
                    help="Number of genomes to add based on top bp intersect.")
parser.add_argument('-e', "--n_extra_genomes_per_species", type=int,
                    default=4,
                    help="Number of genomes to add based on top unique bp intersect.")
parser.add_argument('-t', "--intersect_bp_threshold", type=float,
                    default=500_000,
                    help="Threshold for intersect_bp sourmash result.")
args = parser.parse_args()


def clean_name(s):
    s_ = str(s)
    prefixes = ['GB_', 'GA_', 'RS_']
    suffixes = ['_genomic.fna.gz', '.gff.gz.fna.gz', '.fa.gz', '.fa']
    for v in prefixes:
        if s_.startswith(v):
            s_ = s_.removeprefix(v)
    for v in suffixes:
        if s_.endswith(v):
            s_ = s_.removesuffix(v)
    return s_

genome2species = {}
with open(args.genome_species, 'rt') as f:
    for l in f:
        k,v = [v.strip() for v in l.split('\t')]
        genome2species[k] = v

sourmash_profile = pd.read_csv(args.sourmash_profile)

species_scores = defaultdict(list)
for idx,r in sourmash_profile.iterrows():
    genome = clean_name(r['name'])
    species = genome2species[genome] if genome in genome2species else None
    if r.intersect_bp > args.intersect_bp_threshold:
        species_scores[species].append((r['name'],r.intersect_bp,r.unique_intersect_bp))

genomes = []
for _,vs in species_scores.items():
    add_genomes = {v for v,_,_ in sorted(vs, key=lambda x:-x[1])[:args.n_genomes_per_species]}
    add_genomes |= {v for v,_,_ in sorted(vs, key=lambda x:-x[2])[:args.n_extra_genomes_per_species]}
    genomes.extend(list(add_genomes))

with open(args.output_fp, 'wt') as f:
    for v in genomes:
        f.write(f"{v}\n")