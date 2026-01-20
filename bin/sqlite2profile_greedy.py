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
parser.add_argument('-c', "--genome_cds_filepaths", type=str,
                    required=True,
                    help="File containing list of filepaths for genome to CDS mapping TSV files.")
parser.add_argument('-o', "--output_dir", type=str,
                    required=True,
                    help="Output directory.")
parser.add_argument('-d', "--loc_index_denom", type=int,
                    default=100,
                    help="Denominator for indexing CDS location.")
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


    # load data

    genome2cds_fp = {}
    with open(args.genome_cds_filepaths, 'rt') as f:
        for i,l in enumerate(f):
            k, v = [v.strip() for v in l.split('\t')]
            genome2cds_fp[k] = v


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
        query_idxs = {v for v, in cur}
        
        query_names = set()
        for query_idxs_ in batch(list(query_idxs), query_batch_n):
            cur = db.cursor()
            cur.execute(f'SELECT name FROM query WHERE idx IN ({",".join([str(v) for v in query_idxs_])});')
            query_names.update({f"\"{v}\"" for v, in cur})
        query_idxs = None
        
        reassigned_reads = set()
        for query_names_ in batch(list(query_names), query_batch_n):
            cur = db.cursor()
            cur.execute(f'SELECT idx FROM query WHERE name IN ({",".join(query_names_)});')
            reassigned_reads.update({v for v, in cur})
        query_names = None
        
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
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM species_genome_read_mappings WHERE genome={genome};')
        mappings = defaultdict(lambda :(None,None,None,0,0))
        for q,r,rs,re_,a,ani in cur:
            if a>=mappings[q][-2]:
                mappings[q] = (r,rs,re_,a,ani)
        
        for k,t in mappings.items():
            if t[0] is None:
                continue
            contig_coverage_depth[reference_list[t[0]]][t[1]:t[2]+1] += 1
        
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


    # CDS profile

    cur = db.cursor()
    cur.execute('SELECT species,genome FROM species_genome_read_mappings;')
    species_genomes = defaultdict(set)
    for species,genome in cur:
        species_genomes[species].add(genome)
        
    cluster_list = []
    cluster_index = {}
    species_cds_mappings = {}
    transaction_count = 0
    commit_n = 100_000
    for species, genomes in species_genomes.items():
        
        for genome in genomes:
            if not genome_list[genome] in genome2cds_fp:
                continue
                
            cur = db.cursor()
            cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM species_genome_read_mappings WHERE genome={genome};')
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
                        reference_list.append(contig)
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
                                INSERT INTO cds_read_mappings (
                                    species, cluster, query, partial, genome, reference,
                                    rstart, rend, qstart, qend, cstart, cend, proportion_l
                                )
                                VALUES (
                                    {species}, {cluster}, {k}, {partial}, {genome}, {contig}, 
                                    {t[1]}, {t[2]}, {start_}, {end_}, {start}, {end}, {l/cds_l}
                                );
                            ''')
                            transaction_count += 1

                        if transaction_count%commit_n == 0:
                            db.commit()
                            
            else:
                db.commit()
                        
        # greedy read reassignment
        cur = db.cursor()
        cur.execute(f'SELECT cluster,partial,query,proportion_l FROM cds_read_mappings WHERE species={species};')
        cds_coverages = defaultdict(lambda :defaultdict(int))
        for c,p,q,l in cur:
            if l>cds_coverages[(c,p)][q]:
                cds_coverages[(c,p)][q] = l
        cds_coverages = {k:sum(list(d.values())) for k,d in cds_coverages.items()}
        n_cdss = len({k for k,_ in cds_coverages})
        
        assigned_cdss = set()
        while (len(assigned_cdss)<n_cdss):
            top_cds, top_cds_partial = max(cds_coverages.items(), key=lambda x:(-x[0][1],x[1]))[0]
            assigned_cdss.add(top_cds)
            if (top_cds,0) in cds_coverages:
                del cds_coverages[(top_cds,0)]
            if (top_cds,1) in cds_coverages:
                del cds_coverages[(top_cds,1)]
            # print(datetime.datetime.now(), cluster_list[top_cds], len(assigned_cdss), flush=True)
            
            # get reassigned reads and locations of top cds in all genomes
            cur = db.cursor()
            cur.execute(f'SELECT idx,query,genome,reference,cstart,cend FROM cds_read_mappings WHERE cluster={top_cds};')
            reassigned_reads = set()
            starts_ends = defaultdict(set)
            for idx,q,g,c,cs,ce in cur:
                reassigned_reads.add(q)
                starts_ends[(g,c)].add((cs,ce))
        
            # prefetch reads for reassigned_cdss
            reassigned_cdss = defaultdict(list)
            for reassigned_reads_ in batch(list(reassigned_reads), query_batch_n):
                cur = db.cursor()
                cur.execute(f'''
                    SELECT idx,cluster,query,genome,reference,rstart,rend,qstart,qend 
                    FROM cds_read_mappings 
                    WHERE query IN ({",".join([str(v) for v in reassigned_reads_])});
                ''')
                for idx,cds,q,g,c,rs,re,qs,qe in cur:
                    if cds in assigned_cdss:
                        continue
                    reassigned_cdss[cds].append((idx,q,g,c,rs,re,qs,qe))
                
            for k,vs in reassigned_cdss.items():
                # remove reads that do not overlap with top CDS
                rm_idx = set()
                for idx,q,g,c,rs,re,qs,qe in vs:
                    # start = max([0,qs-rs])
                    # end = min([re-rs,qe-rs])
                    if all([((re<s) or (rs>e)) for s,e in starts_ends[(g,c)]]):
                        rm_idx.add(idx)
                
                # remove reads 
                for rm_idx_ in batch(list(rm_idx), query_batch_n):
                    db.execute(f'DELETE FROM cds_read_mappings WHERE idx IN ({",".join([str(v) for v in rm_idx_])});')
                db.commit()
                
                # recalculate CDS coverage
                cds_coverages[(k,0)] = defaultdict(int)
                cds_coverages[(k,1)] = defaultdict(int)
                cur = db.cursor()
                cur.execute(f'SELECT cluster,partial,query,proportion_l FROM cds_read_mappings WHERE species={species} AND cluster={k};')
                for c,p,q,l in cur:
                    if l>cds_coverages[(c,p)][q]:
                        cds_coverages[(c,p)][q] = l
                cds_coverages = {k_:sum(list(v.values())) if isinstance (v, dict) else v for k_,v in cds_coverages.items()}
        
        print('Assigning reads to CDSs:', datetime.datetime.now(), species_list[species], len(assigned_cdss), flush=True)

    
    species_clusters = defaultdict(set)
    cur = db.cursor()
    cur.execute('SELECT DISTINCT species,cluster FROM cds_read_mappings;')
    for species, cluster in cur:
        species_clusters[species].add(cluster)

    species_cds_cluster_coverage = defaultdict(dict)
    for species, clusters in species_clusters.items():
        print('Calculating CDS coverage:', datetime.datetime.now(), species_list[species], len(clusters), flush=True)
        for cluster in clusters:
            cur = db.cursor()
            cur.execute(f'''
                SELECT query,genome,reference,rstart,rend,cstart,cend,proportion_l 
                FROM cds_read_mappings 
                WHERE (species={species}) AND (cluster={cluster});
            ''')
            mappings = defaultdict(lambda :(None,None,None,None,None,None,0))
            for q,g,r,rs,re,cs,ce,pl in cur:
                if pl>=mappings[q][-1]:
                    mappings[q] = (g,r,rs,re,cs,ce,pl)
            
            cds_coverages = {}
            for _,(g,r,rs,re,cs,ce,_) in mappings.items():
                if not (g,r,cs) in cds_coverages:
                    cds_coverages[(g,r,cs)] = np.zeros(ce-cs)
                cds_coverages[(g,r,cs)][rs-cs:re-cs+1] += 1
                
            cds_coverage_summary = {}
            for k,vs in cds_coverages.items():
                coverage_depth = vs.mean()
                if coverage_depth<=0:
                    continue
                coverage_breadth = (vs>0).mean()
                mean_depth_ = coverage_depth if coverage_depth<700 else 700
                expected_breadth = 1 - (1/(np.log2(1+np.exp(mean_depth_))))
                cds_coverage_summary[k] = (float(coverage_depth), float(coverage_breadth), float(expected_breadth), float(coverage_breadth/expected_breadth), len(vs))

            if cds_coverage_summary:
                species_cds_cluster_coverage[species][cluster] = tuple([max(v) for v in zip(*list(cds_coverage_summary.values()))])
            

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
            query_idxs = {v for v, in cur}

            query_names = set()
            for query_idxs_ in batch(list(query_idxs), query_batch_n):
                cur = db.cursor()
                cur.execute(f'SELECT name FROM query WHERE idx IN ({",".join([str(v) for v in query_idxs_])});')
                query_names.update({f"\"{v}\"" for v, in cur})
            query_idxs = None

            reassigned_reads = set()
            for query_names_ in batch(list(query_names), query_batch_n):
                cur = db.cursor()
                cur.execute(f'SELECT idx FROM query WHERE name IN ({",".join(query_names_)});')
                reassigned_reads.update({v for v, in cur})
            query_names = None

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
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,ani,ani_gapped_fullread FROM species_genome_read_mappings WHERE genome={genome};')
        mappings = defaultdict(set)
        for q,r,rs,re_,a,ani in cur:
            mappings[q].add((r,rs,re_,a,ani))
        
        for k,ts in mappings.items():
            t = max(ts, key=lambda x:x[-2])
            contig_coverage_depth[reference_list[t[0]]][t[1]:t[2]] += 1
        
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
    with gzip.open(out_dir / f"{prefix}species_cds_coverage.tsv.gz", 'wt') as f:
        for species,d in species_cds_cluster_coverage.items():
            for cds,(depth,b,e,r,l) in d.items():
                f.write(f'{species}\t{cluster_list[cds]}\t{depth}\t{b}\t{e}\t{r}\t{l}\n')
    with gzip.open(out_dir / f"{prefix}species_index.txt.gz", 'wt') as f:
        for s in species_list:
            f.write(f'{s}\n')
