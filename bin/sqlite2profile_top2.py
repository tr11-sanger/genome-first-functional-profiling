#!/usr/bin/env python

import numpy as np
import json
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
parser.add_argument('-c', "--genome_cds_filepaths", type=str,
                    required=True,
                    help="File containing list of filepaths for genome to CDS mapping TSV files.")
parser.add_argument('-o', "--output_dir", type=str,
                    required=True,
                    help="Output directory.")
parser.add_argument('-p', "--output_prefix", type=str,
                    default='',
                    help="Output prefix.")
parser.add_argument('-d', "--loc_index_denom", type=int,
                    default=100,
                    help="Denominator for indexing CDS location.")
args = parser.parse_args()

if __name__ == '__main__':

    def batch(l, n):
        for i in range((len(l)//n)+1):
            if i*n<len(l):
                yield l[i*n:(i+1)*n]

    commit_n = 100_000
    query_batch_n = 100_000
    transaction_count = 0


    # load data

    genome2cds_fp = {}
    with open(args.genome_cds_filepaths, 'rt') as f:
        for i,l in enumerate(f):
            k, v = [v.strip() for v in l.split('\t')]
            genome2cds_fp[k] = v


    # load sqlite
    print(datetime.datetime.now(), 'Loading SQLite DB', flush=True)
    
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
    print(datetime.datetime.now(), 'Loading indexes from DB', flush=True)
    
    cur = db.cursor()
    cur.execute(f'''
        SELECT COUNT(DISTINCT ma.query)
        FROM species_genome_read_mappings AS ma
    ''')
    n_mapped_read_ends = int(list(cur)[0][0])
    
    cur = db.cursor()
    cur.execute(f'''
        SELECT COUNT(DISTINCT qu.name)
        FROM species_genome_read_mappings AS ma
        LEFT JOIN query AS qu ON ma.query = qu.idx;
    ''')
    n_mapped_read_pairs = int(list(cur)[0][0])
    
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
    
    
    # Species profile 
    # get top read match per genome
    print(datetime.datetime.now(), 'Trimming mappings to one paired match per genome, and only to genomes of the top matching species', flush=True)
    
    def gen_top_genome_matches():
        cur = db.cursor()
        cur.execute(f'''
            SELECT ma.idx,qu.name,qu.pair,ma.reference,ma.genome,ma.species,ma.ani_gapped_fullread
            FROM species_genome_read_mappings AS ma
            LEFT JOIN query AS qu ON ma.query = qu.idx
            ORDER BY qu.name DESC;
        ''')
    
        def get_good_mappings(paired_matches):
            if not paired_matches:
                return
    
            top_paired_matches = defaultdict(list)
            for (k_gi,k_si),d1 in paired_matches.items():
                top_genome_paired_matches = defaultdict(list)
                for k_ri,d2 in d1.items():
                    if len(d2)==2:
                        for _,l in d2.items():
                            if not l:
                                continue
                            top_genome_paired_matches[k_ri].append(max(l, key=lambda x:x[1]))
                if not top_genome_paired_matches:
                    continue
                top_paired_matches[(k_gi,k_si)] = max(top_genome_paired_matches.items(), key=lambda x:sum([v for _,v in x[1]]))[1]
    
            if not top_paired_matches:
                return
    
            top_species = max(top_paired_matches.items(), key=lambda x:sum([v[1] for v in x[1]]))[0][1]
            for (k_gi,k_si),ms in top_paired_matches.items():
                if k_si != top_species:
                    continue
                for k_mi,_ in ms:
                    yield k_mi
    
        current_qn = None
        paired_matches = defaultdict(lambda :defaultdict(lambda :defaultdict(list)))
        for mi,qn,qp,ri,gi,si,ani in cur:
            if qn != current_qn:
                for v in get_good_mappings(paired_matches):
                    yield v
                paired_matches = defaultdict(lambda :defaultdict(lambda :defaultdict(list)))
            paired_matches[(gi,si)][ri][qp].append((mi,ani))
            current_qn = qn
        else:
            for v in get_good_mappings(paired_matches):
                yield v
    
    good_mappings = []
    for mi in gen_top_genome_matches():
        good_mappings.append(mi)
        if len(good_mappings)>=query_batch_n:
            db.execute(f'''
                INSERT INTO top_species_genome_read_mappings (species, genome, query, reference, rstart, rend, ani, ani_gapped, ani_gapped_fullread)
                SELECT species, genome, query, reference, rstart, rend, ani, ani_gapped, ani_gapped_fullread 
                FROM species_genome_read_mappings
                WHERE idx IN ({",".join([str(v) for v in good_mappings])});
            ''')
            transaction_count += len(good_mappings)
            good_mappings = []
    
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
        print(datetime.datetime.now(), 'Calculating genome coverage', species_list[species], genome_list[genome], flush=True)
    
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
    
        mean_ani = np.mean([v[4] for _,v in mappings.items()])
    
        species_genomes_coverage[species][genome] = (
            float(genome_coverage_depth), 
            float(genome_coverage_breadth), 
            float(genome_expected_breadth), 
            float(genome_coverage_breadth/genome_expected_breadth), 
            float(mean_ani), 
            len(mappings), 
            mapped_read_pairs
        )
    
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
    mapped_reads = {k:len(v) for k,v in mapped_reads.items()}
    mapped_read_ends = {k:len(v) for k,v in mapped_read_ends.items()}
    
    species_top_genome_coverage = {k:sorted(d.items(), key=lambda x:x[1][1])[-1][1] for k,d in species_genomes_coverage.items()}
    species_top_genome_coverage = {k:v[:-2]+(mapped_read_ends[k],mapped_reads[k]) for k,v in species_top_genome_coverage.items()}
    
    
    # genome coverage
    genomes_coverage = {k:v for _,d in species_genomes_coverage.items() for k,v in d.items()}
    

    # CDS profile
    species_cds_cluster_coverage = defaultdict(dict)

    cur = db.cursor()
    cur.execute('SELECT DISTINCT species,genome FROM top_species_genome_read_mappings;')
    species_genomes = defaultdict(set)
    for species,genome in cur:
        species_genomes[species].add(genome)
        
    cluster_list = []
    cluster_index = {}
    species_cds_mappings = {}
    for species, genomes in species_genomes.items():
        
        for genome in genomes:
            if not genome_list[genome] in genome2cds_fp:
                continue
                
            print(datetime.datetime.now(), 'Mapping to CDSs', species_list[species], genome_list[genome], flush=True)

            cur = db.cursor()
            cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM top_species_genome_read_mappings WHERE genome={genome};')
            mappings = defaultdict(set)
            for q,r,rs,re_,a,ani in cur:
                mappings[q].add((r,rs,re_,a,ani))

            cds_file = []
            with open(genome2cds_fp[genome_list[genome]], 'rt') as f:
                for l in f:
                    name, cluster, contig, _, start, end, direction, data = [v.strip() for v in l.split('\t')]
                    
                    if cluster not in cluster_index:
                        cluster_index[cluster] = len(cluster_list)
                        cluster_list.append(cluster)
                        db.execute(f'''
                            INSERT INTO cluster (name)
                            VALUES ("{cluster}");
                        ''')
                        transaction_count += 1
                        
                        if transaction_count%commit_n == 0:
                            db.commit()
                    
                    if contig not in reference_index:
                        reference_index[contig] = len(reference_list)
                        reference_list[len(reference_list)] = contig
                        db.execute(f'''
                            INSERT INTO reference (name)
                            VALUES ("{contig}");
                        ''')
                        transaction_count += 1
                        
                        if transaction_count%commit_n == 0:
                            db.commit()
                    
                    data_dict = {}
                    for s in data.split(';'):
                        k,v = s.split('=')
                        data_dict[k] = v
                    partial = data_dict['partial']!='00'
                    cds_file.append((cluster_index[cluster], reference_index[contig], int(start), int(end), partial))
                    
            cds_file_index = defaultdict(lambda :defaultdict(set))
            for i,(cluster, contig, start, end, partial) in enumerate(cds_file):
                for j in range(start//args.loc_index_denom, (end+1)//args.loc_index_denom):
                    cds_file_index[contig][j].add(i)

            # translate mappings to CDSs
            for k,ts in mappings.items():
                for t in ts:
                    if t[0] not in cds_file_index:
                        continue
                    
                    cdss = {v for i in range(t[1]//args.loc_index_denom,(t[2]+1)//args.loc_index_denom) for v in cds_file_index[t[0]][i]}
                    
                    for cds_i in cdss:
                        cluster, contig, start, end, partial = cds_file[cds_i]
                        start_ = max([t[1], start])
                        end_ = min([t[2], end])
                        l = end_-start_
                        cds_l = end-start
                        if l>0:
                            db.execute(f'''
                                INSERT INTO species_cds_read_mappings (
                                    species, cluster, query, partial, genome, reference,
                                    rstart, rend, qstart, qend, cstart, cend, proportion_l, ani_gapped_fullread
                                )
                                VALUES (
                                    {species}, {cluster}, {k}, {partial}, {genome}, {contig}, 
                                    {t[1]}, {t[2]}, {start_}, {end_}, {start}, {end}, {l/cds_l}, {t[4]}
                                );
                            ''')
                            transaction_count += 1

                        if transaction_count%commit_n == 0:
                            db.commit()
                            
            else:
                db.commit()
                        
        # top read assignment
        cur = db.cursor()
        cur.execute(f'SELECT cluster,partial,query,genome,reference,proportion_l,rstart,rend,cstart,cend FROM species_cds_read_mappings ORDER BY query, ani_gapped_fullread DESC;')
        current_read = None
        read_covered_bases = None
        cds_coverages = {}
        for c,p,q,g,r,l,rs,re_,cs,ce in cur:
            if not q==current_read:
                read_covered_bases = set()
                current_read = q
            
            read_mapping_bases = set(range(rs-cs, re_-cs+1))
            if read_mapping_bases & read_covered_bases:
                continue
            read_covered_bases.update(read_mapping_bases)

            if not (g,r,c,cs) in cds_coverages:
                cds_coverages[(g,r,c,cs)] = np.zeros(ce-cs)
            cds_coverages[(g,r,c,cs)][rs-cs:re_-cs+1] += 1

        cds_coverage_summary = defaultdict(list)
        for k,v in cds_coverages.items():
            coverage_depth = v.mean()
            if coverage_depth<=0:
                continue
            coverage_breadth = (v>0).mean()
            mean_depth_ = coverage_depth if coverage_depth<700 else 700
            expected_breadth = 1 - (1/(np.log2(1+np.exp(mean_depth_))))
            cds_coverage_summary[k[2]].append((float(coverage_depth), float(coverage_breadth), float(expected_breadth), float(coverage_breadth/expected_breadth), len(v)))

        for cluster, l in cds_coverage_summary.items():
            ds, bs, es, rs, ls = zip(*list(l))
            b, e, l_ = max(zip(bs, es, ls), key=lambda x:x[0])
            species_cds_cluster_coverage[species][cluster] = (sum(ds), b, e, b/e, l_)

        db.execute(f'''
            DELETE FROM species_cds_read_mappings; 
        ''')
        db.commit()
    

    # mapping statistics
    print(datetime.datetime.now(), f'Generating mapping statistics', flush=True)
    
    cur = db.cursor()
    cur.execute(f'''
        SELECT COUNT(idx)
        FROM query;
    ''')
    n_read_ends = int(list(cur)[0][0])
    
    cur = db.cursor()
    cur.execute(f'''
        SELECT COUNT(DISTINCT name)
        FROM query;
    ''')
    n_read_pairs = int(list(cur)[0][0])
    
    cur = db.cursor()
    cur.execute(f'''
        SELECT COUNT(DISTINCT qu.name)
        FROM top_species_genome_read_mappings AS ma
        LEFT JOIN query AS qu ON ma.query = qu.idx;
    ''')
    n_paired_mapped_reads = int(list(cur)[0][0])
    
    mapping_statistics = {
        'n_read_ends': n_read_ends,
        'n_read_pairs': n_read_pairs,
        'n_mapped_read_ends': n_mapped_read_ends,
        'n_mapped_read_pairs': n_mapped_read_pairs,
        'n_paired_mapped_reads': n_paired_mapped_reads,
    }
    
    
    # Outputs
    print(datetime.datetime.now(), f'Writing outputs', flush=True)
    
    prefix = f"{args.output_prefix}_" if len(args.output_prefix)>0 else ""
    out_dir = Path(args.output_dir)
    os.makedirs(out_dir, exist_ok=True)
    
    with gzip.open(out_dir / f"{prefix}species_coverage.tsv.gz", 'wt') as f:
        for species,(d,b,e,r,a,n1,n2) in species_top_genome_coverage.items():
            f.write(f'{species_list[species]}\t{d}\t{b}\t{e}\t{r}\t{a}\t{n1}\t{n2}\n')

    with gzip.open(out_dir / f"{prefix}genome_coverage.tsv.gz", 'wt') as f:
        for genome,(d,b,e,r,a,n1,n2) in genomes_coverage.items():
            f.write(f'{genome_list[genome]}\t{d}\t{b}\t{e}\t{r}\t{a}\t{n1}\t{n2}\n')

    with gzip.open(out_dir / f"{prefix}species_cds_coverage.tsv.gz", 'wt') as f:
        for species,d in species_cds_cluster_coverage.items():
            for cds,(depth,b,e,r,l) in d.items():
                f.write(f'{species}\t{cluster_list[cds]}\t{depth}\t{b}\t{e}\t{r}\t{l}\n')

    with gzip.open(out_dir / f"{prefix}species_index.tsv.gz", 'wt') as f:
        for i,s in species_list.items():
            f.write(f'{i}\t{s}\n')

    with open(out_dir / f"{prefix}mapping_statistics.json", 'wt') as f:
        json.dump(mapping_statistics, f)
    
    print(datetime.datetime.now(), f'Complete', flush=True)
