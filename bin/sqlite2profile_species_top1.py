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

parser = argparse.ArgumentParser(description='Read SQLite database containing read mapping information to generate taxonomic and functional profiles. Calculates coverage depth and breadth for each species, genome and CDS.')
parser.add_argument('-s', "--sqlite", type=str,
                    required=True,
                    help="SQLite database path (GZiped) containing read mapping data.")
parser.add_argument('-o', "--output_dir", type=str,
                    required=True,
                    help="Output directory.")
parser.add_argument('-p', "--output_prefix", type=str,
                    default='',
                    help="Output prefix.")
parser.add_argument('-t', "--min_coverage_ratio", type=float,
                    default=0.95,
                    help="Minimum observed:expected coverage breadth ratio to exclude species and reassign their reads.")
args = parser.parse_args()

if __name__ == '__main__':

    def batch(l, n):
        for i in range((len(l)//n)+1):
            if i*n<len(l):
                yield l[i*n:(i+1)*n]

    commit_n = 100_000
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
        SELECT COUNT(idx)
        FROM species_genome_read_mappings;
    ''')
    n_queries = int(list(cur)[0][0])

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
    print('Assigning reads to species:', datetime.datetime.now(), flush=True)
    all_assigned = False
    while not all_assigned:
        # get top read match per genome
        print('Trimming mappings to top match per genome:', datetime.datetime.now(), flush=True)

        def gen_top_matches():
            cur = db.cursor()
            cur.execute(f'''
                SELECT idx,query,genome,species,ani_gapped_fullread
                FROM species_genome_read_mappings
                ORDER BY query,genome,ani_gapped_fullread DESC;
            ''')

            current_qi = None
            current_gi = None
            for mi,qi,gi,si,ani in cur:
                if (gi != current_gi) or (qi != current_qi):
                    yield mi,qi,si
                    
                current_gi = gi
                current_qi = qi

        # identify mappings that map to the top species
        print('Trimming mappings to only match to top species:', datetime.datetime.now(), flush=True)
        assigned_species = np.ones(n_queries, dtype=int) * -1
        good_mappings = []
        for i,q,s in gen_top_matches():
            if assigned_species[q]==-1:
                assigned_species[q] = s
                good_mappings.append(i)
            else:
                if s == assigned_species[q]:
                    good_mappings.append(i)

        # populate top-mappings table
        for good_mappings_  in batch(list(good_mappings), query_batch_n):
            db.execute(f'''
                INSERT INTO top_species_genome_read_mappings (species, genome, query, reference, rstart, rend, ani, ani_gapped, ani_gapped_fullread)
                SELECT species, genome, query, reference, rstart, rend, ani, ani_gapped, ani_gapped_fullread 
                FROM species_genome_read_mappings
                WHERE idx IN ({",".join([str(v) for v in good_mappings_])});
            ''')
            transaction_count += len(good_mappings_)

            if transaction_count>=commit_n:
                transaction_count = 0
                db.commit()
        else:
            db.commit()
    
        # calculate genome coverages
        species_genomes_coverage = defaultdict(dict)
        cur = db.cursor()
        cur.execute('SELECT DISTINCT species,genome FROM top_species_genome_read_mappings;')
        species_genome = [(s,g) for s,g in cur]
        for species,genome in species_genome:
            print('Calculating species coverage:', datetime.datetime.now(), species_list[species], genome_list[genome], flush=True)
            
            contig_coverage_depth = {v: np.zeros(contig_lengths[v], dtype=int) for v in genome2contigs[genome]}
            
            cur = db.cursor()
            cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM top_species_genome_read_mappings WHERE genome={genome};')
            mappings = {}
            for q,r,rs,re_,a,ani in cur:
                mappings[q] = (r,rs,re_,a,ani)
            
            for k,t in mappings.items():
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
        mapped_read_pairs = None
        contig_coverage_depth = None

        # species read counts
        cur = db.cursor()
        cur.execute(f'SELECT ma.species, qu.idx, qu.name FROM top_species_genome_read_mappings as ma LEFT JOIN query as qu ON ma.query=qu.idx;')
        mapped_reads = defaultdict(set)
        mapped_read_ends = defaultdict(set)
        for s,qi,qn in cur:
            mapped_read_ends[s].add(qi)
            mapped_reads[s].add(qn)
         
        species_top_genome_coverage = {k:sorted(d.items(), key=lambda x:x[1][1])[-1][1] for k,d in species_genomes_coverage.items()}
        species_top_genome_coverage = {k:v[:-2]+(len(mapped_read_ends[k]),len(mapped_reads[k])) for k,v in species_top_genome_coverage.items()}

        # remove bottom species, or all species that fail
        present_species = {str(k) for k,(d,b,e,r,n,n_) in species_top_genome_coverage.items() if r>=args.min_coverage_ratio}
    
        if len(present_species) < len(species_top_genome_coverage):
            print(f'Removing {len(species_top_genome_coverage) - len(present_species)} species that are below coverage threshold:', datetime.datetime.now(), flush=True)
            # delete mappings not in present species
            db.execute(f'DELETE FROM top_species_genome_read_mappings;')
            db.execute(f'DELETE FROM species_genome_read_mappings WHERE species NOT IN ({",".join(present_species)});')
        else:
            all_assigned = True
            # remove non-top read mapping data
            db.execute(f'DELETE FROM species_genome_read_mappings;')
        db.commit()

    # genome coverage
    genomes_coverage = {k:v for _,d in species_genomes_coverage.items() for k,v in d.items()}

    print(f'Writing outputs', datetime.datetime.now(), flush=True)

    # Outputs
    prefix = f"{args.output_prefix}_" if len(args.output_prefix)>0 else ""
    out_dir = Path(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)

    with gzip.open(out_dir / f"{prefix}species_coverage.tsv.gz", 'wt') as f:
        for species,(d,b,e,r,n1,n2) in species_top_genome_coverage.items():
            f.write(f'{species_list[species]}\t{d}\t{b}\t{e}\t{r}\t{n1}\t{n2}\n')
    with gzip.open(out_dir / f"{prefix}genome_coverage.tsv.gz", 'wt') as f:
        for genome,(d,b,e,r,n1,n2) in genomes_coverage.items():
            f.write(f'{genome_list[genome]}\t{d}\t{b}\t{e}\t{r}\t{n1}\t{n2}\n')

    print(f'Complete', datetime.datetime.now(), flush=True)