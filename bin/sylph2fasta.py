import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(description='Choose reference genomes based on Sylph profile.')
parser.add_argument('-s', "--sylph_profile", type=str,
                    required=True,
                    help="Sylph profile filepath.")
parser.add_argument('-o', "--output_fp", type=str,
                    required=True,
                    help="Output filepath.")
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


sylph_profile = pd.read_csv(args.sylph_profile, sep='\t')

genomes = [clean_name(v.split('/')[-1]) for v in sylph_profile.Genome_file]

with open(args.output_fp, 'wt') as f:
    for v in genomes:
        f.write(f"{v}\n")
