#!/usr/bin/env python

import numpy as np
import re
import gzip
from pathlib import Path
from collections import defaultdict
import datetime
import sys
import os
import argparse
import sqlite3
import subprocess

parser = argparse.ArgumentParser(description='Read SQLite database containing read mapping information to generate taxonomic and functional profiles. Reassigns multiple mappings with winner-takes-all strategy. Calculates coverage depth and breadth for each species, genome and CDS.')
parser.add_argument('-s', "--sqlite", type=str,
                    required=True,
                    help="SQLite database path (GZiped) containing read mapping data.")
parser.add_argument('-o', "--output_dir", type=str,
                    required=True,
                    help="Output directory.")
parser.add_argument('-p', "--output_prefix", type=str,
                    default='',
                    help="Output prefix.")
args = parser.parse_args()


if __name__ == '__main__':

    def batch(l, n):
        for i in range((len(l)//n)+1):
            if i*n<len(l):
                yield l[i*n:(i+1)*n]

    query_batch_n = 100_000
    transaction_count = 0


    # load sqlite

    if args.sqlite[-3:]=='.gz':
        backup_db_path = args.sqlite[:-3]
        subprocess.call(f"gunzip -c {args.sqlite} > {backup_db_path}", shell=True)
    else:
        backup_db_path = args.sqlite

    db = sqlite3.connect(':memory:')
    backup_db = sqlite3.connect(backup_db_path)
    with backup_db:
        backup_db.backup(db)
    backup_db.close()

    if args.sqlite[-3:]=='.gz':
        os.remove(backup_db_path)


    # load indexes

    cur = db.cursor()
    cur.execute(f'''
        SELECT idx, name
        FROM species;
    ''')
    species_list = {idx: species_name for idx, species_name in cur}

    cur = db.cursor()
    cur.execute(f'''
        SELECT idx, name
        FROM genome;
    ''')
    genome_list = {idx: genome_name for idx, genome_name in cur}

    cur = db.cursor()
    cur.execute(f'''
        SELECT idx, name, genome, length
        FROM reference;
    ''')
    contig_lengths = {}
    genome2contigs = defaultdict(set)
    reference_list = {}
    reference_index = {}
    for idx, contig_name, genome_idx, length in cur:
        contig_lengths[idx] = length
        genome2contigs[genome_idx].add(idx)
        reference_list[idx] = contig_name
        reference_index[contig_name] = idx
    genome2contigs = dict(genome2contigs)

 
    # filter reads by pairing
    print('Filtering unpaired read mappings:', datetime.datetime.now(), flush=True)
    cur = db.cursor()
    cur.execute(f'''
        SELECT qu.name, qu.pair, ma.idx, ma.ani_gapped_fullread, ma.reference, ma.genome, ma.species
        FROM species_genome_read_mappings AS ma
        LEFT JOIN query AS qu ON ma.query = qu.idx
        ORDER BY qu.name;'''
    )
    
    current_qn = None
    paired_mappings = defaultdict(lambda :defaultdict(set))
    unpaired_mappings = []
    
    def remove_unpaired_from_db(qn, paired_mappings, force_commit=True):
        global unpaired_mappings
        global transaction_count
        
        a = set(paired_mappings[0].keys())
        b = set(paired_mappings[1].keys())
        unpaired_refs = (a-b) | (b-a)
        unpaired_mappings_ = list({m for _,d in paired_mappings.items() for k in unpaired_refs for m in d[k]})
        unpaired_mappings += unpaired_mappings_
        transaction_count += len(unpaired_mappings_)
        
        if force_commit or (transaction_count>query_batch_n):
            db.execute(f'DELETE FROM species_genome_read_mappings WHERE idx IN ({",".join([str(v) for v in unpaired_mappings])});')
            db.commit()
            unpaired_mappings = []
            transaction_count = 0
    
    for qn,qp,m,ani,r,g,s in cur:
        if current_qn is None:
            current_qn = qn
        if current_qn != qn:
            remove_unpaired_from_db(current_qn, paired_mappings, force_commit=False)
            paired_mappings = defaultdict(lambda :defaultdict(set))
            current_qn = qn
        paired_mappings[qp][s].add(m)
    else:
        remove_unpaired_from_db(current_qn, paired_mappings, force_commit=True)
        del paired_mappings


    # Species profile 

    assigned_species = set()

    print('Assigning reads to species:', datetime.datetime.now(), len(assigned_species), flush=True)
    
    cur = db.cursor()
    cur.execute('SELECT COUNT(DISTINCT species) FROM species_genome_read_mappings;')
    n_species = int(cur.fetchone()[0])

    while len(assigned_species)<n_species:
        cur = db.cursor()
        cur.execute(f'SELECT species,genome,query FROM species_genome_read_mappings WHERE species NOT IN ({",".join([str(v) for v in assigned_species])});')
        species_read_counts = defaultdict(lambda :defaultdict(set))
        for s,g,q in cur:
            species_read_counts[s][g].add(q)
        species_read_counts = {k:max([len(vs) for vs in d.values()]) for k,d in species_read_counts.items()}
        
        if len(species_read_counts)==0:
            break

        top_species = max(species_read_counts.items(), key=lambda x:x[1])[0]
        assigned_species.add(top_species)
        del species_read_counts[top_species]
        
        print('Assigning reads to species:', datetime.datetime.now(), len(assigned_species), n_species, species_list[top_species], flush=True)
        
        cur = db.cursor()
        cur.execute(f'SELECT query FROM species_genome_read_mappings WHERE species={top_species};')
        reassigned_reads = {v for v, in cur}
        
        for reassigned_reads_ in batch(list(reassigned_reads), query_batch_n):
            db.execute(f'DELETE FROM species_genome_read_mappings WHERE (query IN ({",".join([str(v) for v in reassigned_reads_])})) AND (species!={top_species});')
        reassigned_reads = None
        db.commit()    
        
    species_read_counts = None

    species_genomes_coverage = defaultdict(dict)
    cur = db.cursor()
    cur.execute('SELECT DISTINCT species,genome FROM species_genome_read_mappings;')
    for species, genome in cur:
        print('Calculating species coverage:', datetime.datetime.now(), species_list[species], flush=True)
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v], dtype=int) for v in genome2contigs[genome]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM species_genome_read_mappings WHERE genome={genome};')
        mappings = defaultdict(lambda :(None,None,None,0,0))
        for q,r,rs,re_,a,ani in cur:
            if ani>=mappings[q][-1]:
                mappings[q] = (r,rs,re_,a,ani)
        
        for k,t in mappings.items():
            if t[0] is None:
                continue
            contig_coverage_depth[t[0]][t[1]:t[2]+1] += 1
        
        sum_len = sum([len(v) for v in contig_coverage_depth.values()])
        if sum_len==0:
            continue
        genome_coverage_depth = sum([v.sum() for v in contig_coverage_depth.values()]) / sum_len
        genome_coverage_breadth = sum([(v>0).sum() for v in contig_coverage_depth.values()]) / sum_len

        mean_depth_ = genome_coverage_depth if genome_coverage_depth<700 else 700
        genome_expected_breadth = 1 - (1/(np.log2(1+np.exp(mean_depth_))))  # * np.log(1+np.exp(0))))

        mapped_read_pairs = set()
        for mappings_ in batch(list(mappings), query_batch_n):
            cur = db.cursor()
            cur.execute(f'SELECT name FROM query WHERE idx IN ({",".join([str(v) for v in mappings_])});')
            mapped_read_pairs.update({v for v, in cur})
        mapped_read_pairs = len(mapped_read_pairs)

        species_genomes_coverage[species][genome] = (float(genome_coverage_depth), float(genome_coverage_breadth), float(genome_expected_breadth), float(genome_coverage_breadth/genome_expected_breadth), len(mappings), mapped_read_pairs)
    
    mappings = None
    contig_coverage_depth = None
        
    species_top_genome_coverage = {k:sorted(d.items(), key=lambda x:x[1][1])[-1] for k,d in species_genomes_coverage.items()}


    # Genome profile

    cur = db.cursor()
    cur.execute('SELECT species,genome FROM species_genome_read_mappings;')
    species_genomes = defaultdict(set)
    for species,genome in cur:
        species_genomes[species].add(genome)

    for species, genomes in species_genomes.items():
        assigned_genomes = set()
        
        while (len(assigned_genomes)<len(genomes)):
            cur = db.cursor()
            cur.execute(f'SELECT genome,query FROM species_genome_read_mappings WHERE (species={species}) AND (genome NOT IN ({",".join([str(v) for v in assigned_genomes])}));')
            genome_read_counts = defaultdict(set)
            for k,v in cur:
                genome_read_counts[k].add(v)
            genome_read_counts = {k:len(v) for k,v in genome_read_counts.items()}

            if len(genome_read_counts)==0:
                break
            
            top_genome = max(genome_read_counts.items(), key=lambda x:x[1])[0]
            assigned_genomes.add(top_genome)
            del genome_read_counts[top_genome]
            
            print('Assigning reads to genomes:', datetime.datetime.now(), len(assigned_genomes), species_list[species], genome_list[top_genome], flush=True)
            
            cur = db.cursor()
            cur.execute(f'SELECT query FROM species_genome_read_mappings WHERE genome={top_genome};')
            reassigned_reads = {v for v, in cur}

            for reassigned_reads_ in batch(list(reassigned_reads), query_batch_n):
                db.execute(f'DELETE FROM species_genome_read_mappings WHERE (query IN ({",".join([str(v) for v in reassigned_reads_])})) AND (genome!={top_genome});')
            reassigned_reads = None
            db.commit()

    genome_read_counts = None 

    genomes_coverage = {}
    cur = db.cursor()
    cur.execute('SELECT DISTINCT species,genome FROM species_genome_read_mappings;')
    species_genome = [(s,g) for s,g in cur]
    for species,genome in species_genome:
        print('Calculating genome coverage:', datetime.datetime.now(), species_list[species], genome_list[genome], flush=True)
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v], dtype=int) for v in genome2contigs[genome]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM species_genome_read_mappings WHERE genome={genome};')
        mappings = defaultdict(set)
        for q,r,rs,re_,a,ani in cur:
            mappings[q].add((r,rs,re_,a,ani))
        
        for k,ts in mappings.items():
            t = max(ts, key=lambda x:x[-1])
            contig_coverage_depth[t[0]][t[1]:t[2]] += 1
        
        sum_len = sum([len(v) for v in contig_coverage_depth.values()])
        if sum_len==0:
            continue
        genome_coverage_depth = sum([v.sum() for v in contig_coverage_depth.values()]) / sum_len
        genome_coverage_breadth = sum([(v>0).sum() for v in contig_coverage_depth.values()]) / sum_len

        mean_depth_ = genome_coverage_depth if genome_coverage_depth<700 else 700
        genome_expected_breadth = 1 - (1/(np.log2(1+np.exp(mean_depth_))))  # * np.log(1+np.exp(0))))

        mapped_read_pairs = set()
        for mappings_ in batch(list(mappings), query_batch_n):
            cur = db.cursor()
            cur.execute(f'SELECT name FROM query WHERE idx IN ({",".join([str(v) for v in mappings_])});')
            mapped_read_pairs.update({v for v, in cur})
        mapped_read_pairs = len(mapped_read_pairs)

        genomes_coverage[genome] = (float(genome_coverage_depth), float(genome_coverage_breadth), float(genome_expected_breadth), float(genome_coverage_breadth/genome_expected_breadth), len(mappings), mapped_read_pairs)
    
    mappings = None


    # Outputs
    prefix = f"{args.output_prefix}_" if len(args.output_prefix)>0 else ""
    out_dir = Path(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    with gzip.open(out_dir / f"{prefix}species_coverage.tsv.gz", 'wt') as f:
        for species,(genome,(d,b,e,r,n1,n2)) in species_top_genome_coverage.items():
            f.write(f'{species_list[species]}\t{d}\t{b}\t{e}\t{r}\t{n1}\t{n2}\n')
    with gzip.open(out_dir / f"{prefix}genome_coverage.tsv.gz", 'wt') as f:
        for genome,(d,b,e,r,n1,n2) in genomes_coverage.items():
            f.write(f'{genome_list[genome]}\t{d}\t{b}\t{e}\t{r}\t{n1}\t{n2}\n')
