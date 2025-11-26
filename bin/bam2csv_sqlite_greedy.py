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

parser = argparse.ArgumentParser(description='Parse BAM to generate taxonomic and functional profiles. Reassigns multiple mappings with winner-takes-all strategy. Calculates coverage depth and breadth for each species, genome and CDS.')
parser.add_argument('-r', "--refs", type=str,
                    required=True,
                    help="Reference genomes FASTA filepath.")
parser.add_argument('-g', "--genome_species", type=str,
                    required=True,
                    help="Genome to species mapping TSV filepath.")
parser.add_argument('-m', "--genome_contig_mapping", type=str,
                    required=True,
                    help="Genome to contig mapping TSV filepath.")
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
parser.add_argument('-k', "--min_align_length", type=int,
                    default=142,
                    help="Minimum alignment length.")
parser.add_argument('-l', "--max_align_length", type=int,
                    default=158,
                    help="Maximum alignment length.")
parser.add_argument('-a', "--min_ani", type=float,
                    default=0.95,
                    help="Minimum ANI for read mapping.")
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

def read_fasta(fp):
    seqs = {}
    header = None
    seq_lines = []
    with open(fp, 'rt') as f:
        for l in f:
            if l[0]=='>':
                if not header is None:
                    seqs[header] = ''.join(seq_lines)
                
                header = l.strip()[1:].split()[0]
                seq_lines = []
            else:
                if l.strip():
                    seq_lines.append(l.strip())
        else:
            seqs[header] = ''.join(seq_lines)
    return seqs

def batch(l, n):
    for i in range((len(l)//n)+1):
        if i*n<len(l):
            yield l[i*n:(i+1)*n]


if __name__ == '__main__':
    # load data

    genome2cds_fp = {}
    with open(args.genome_cds_filepaths, 'rt') as f:
        for i,l in enumerate(f):
            k, v = [v.strip() for v in l.split('\t')]
            genome2cds_fp[k] = v

    genome2contigs = defaultdict(set)
    with open(args.genome_contig_mapping, 'rt') as f:
        for l in f:
            a,b = [v.strip() for v in l.split(',')]
            genome2contigs[clean_name(a)].add(b)
    contig2genome = {v:k for k,vs in genome2contigs.items() for v in vs}
    genome2contigs = dict(genome2contigs)

    genome2species = {}
    with open(args.genome_species, 'rt') as f:
        for l in f:
            k,v = [v.strip() for v in l.split('\t')]
            genome2species[k] = v

    contig_lengths = {k:len(v) for k,v in read_fasta(args.refs).items()}


    # Set up database

    db_path = ':memory:'  # 
    db = sqlite3.connect(db_path)

    db.execute('''
        CREATE TABLE species (
            idx INTEGER PRIMARY KEY,
            name TEXT
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE genome (
            idx INTEGER PRIMARY KEY,
            name TEXT
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE query (
            idx INTEGER PRIMARY KEY,
            name TEXT,
            pair INTEGER
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE reference (
            idx INTEGER PRIMARY KEY,
            name TEXT
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE cluster (
            idx INTEGER PRIMARY KEY,
            name TEXT
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE species_genome_read_mappings (
            idx INTEGER PRIMARY KEY,
            species INTEGER,
            genome INTEGER,
            query INTEGER,
            reference INTEGER,
            rstart INTEGER,
            rend INTEGER, 
            align_score REAL, 
            ani REAL
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE cds_read_mappings (
            idx INTEGER PRIMARY KEY,
            species INTEGER,
            cluster INTEGER,
            query INTEGER,
            partial INTEGER,
            genome INTEGER,
            reference INTEGER,
            rstart INTEGER,
            rend INTEGER, 
            qstart INTEGER,
            qend INTEGER,
            cstart INTEGER,
            cend INTEGER,
            proportion_l REAL
        );
    ''')
    db.commit()


    # Read BAM

    cigar_codes = {
        'M': (True, True),
        'I': (True, False),
        'D': (False, True),
        'N': (False, True),
        'S': (True, False),
        'H': (False, False),
        'P': (False, False),
        '=': (True, True),
        'X': (True, True),
    }

    def parse_cigar(cigar_str):
        cigar_codes = {
            'M': (True, True),
            'I': (True, False),
            'D': (False, True),
            'N': (False, True),
            'S': (True, False),
            'H': (False, False),
            'P': (False, False),
            '=': (True, True),
            'X': (True, True),
        }
        qalength = 0
        ralength = 0
        alength = 0
        consumes = (False, False)
        v = []
        for c in cigar_str:
            if c in cigar_codes:
                v_int = int(''.join(v))
                consumes = cigar_codes[c]
                if consumes[0]:
                    qalength += v_int
                if consumes[1]:
                    ralength += v_int
                if consumes[0] or consumes[1]:
                    alength += v_int
                v = []
            else:
                v.append(c)
    
        return qalength, ralength, alength
        
    query_count = 0
    query_index = {}
    reference_list = []
    reference_index = {}
    genome_list = []
    genome_index = {}
    species_list = []
    species_index = {}

    transaction_count = 0

    commit_n = 100_000
    read_batch_n = 10_000_000
    query_batch_n = 100_000
    nm_suffix = 'NM:i:'
    as_suffix = 'AS:i:'

    def parse_sam_line(read):
        global query_count
        global query_index
        global reference_list
        global reference_index
        global genome_list
        global genome_index
        global species_list
        global species_index
        global transaction_count
        global commit_n
        global read_batch_n
        global query_batch_n
        global nm_suffix
        global as_suffix
        global contig2genome
        global genome2species

        read_ = read.split('\t')
        
        if not read_[9]:
            return
    
        for v in read_:
            if v[:len(nm_suffix)]==nm_suffix:
                nm_tag = int(v[len(nm_suffix):])
                continue
            if v[:len(as_suffix)]==as_suffix:
                as_tag = int(v[len(as_suffix):])
                continue
    
        mapq = int(read_[4])
        if mapq==0:
            return
            
        query = read_[0]
        _, ralength, alength = parse_cigar(read_[5])
        mapq = int(read_[4])
        flags_int = int(read_[1])
        paired = flags_int & 1
        read1 = flags_int & 64
        read2 = (False if read1 else True) if paired else None
        # forward = flags_int & (16 if read1 else 32) 
        reference = read_[2]
        # rlength = int(read_[7])
        rstart = int(read_[3])
        rend = rstart+ralength
        ani = 1-(nm_tag/alength)
        align_score = as_tag

        if not all([ani>=args.min_ani, alength>=args.min_align_length, alength<=args.max_align_length]):
            return

        genome = contig2genome[reference]
        species = genome2species[genome] if genome in genome2species else 'unknown'

        if not (query,read2) in query_index:
            query_index[(query,read2)] = int(query_count)
            query_count += 1
            db.execute(f'''
                INSERT INTO query (name, pair) 
                VALUES ("{query}", {read2});
            ''')
            transaction_count += 1
            if transaction_count%commit_n == 0:
                db.commit()
        
        if not reference in reference_index:
            reference_index[reference] = len(reference_list)
            reference_list.append(reference)
            db.execute(f'''
                INSERT INTO reference (name) 
                VALUES ("{reference}");
            ''')
            transaction_count += 1
            if transaction_count%commit_n == 0:
                db.commit()
        
        if not genome in genome_index:
            genome_index[genome] = len(genome_list)
            genome_list.append(genome)
            db.execute(f'''
                INSERT INTO genome (name) 
                VALUES ("{genome}");
            ''')
            transaction_count += 1
            if transaction_count%commit_n == 0:
                db.commit()
        
        if not species in species_index:
            species_index[species] = len(species_list)
            species_list.append(species)
            db.execute(f'''
                INSERT INTO species (name)
                VALUES ("{species}");
            ''')
            transaction_count += 1
            if transaction_count%commit_n == 0:
                db.commit()
            
        db.execute(f'''
            INSERT INTO species_genome_read_mappings (
                species, genome, query, reference,
                rstart, rend, align_score, ani
            ) 
            VALUES (
                {species_index[species]}, {genome_index[genome]}, {query_index[(query,read2)]}, {reference_index[reference]}, 
                {rstart}, {rend}, {align_score}, {ani}
            );
        ''')
        transaction_count += 1
        
        if transaction_count%commit_n == 0:
            db.commit()

    reads = []
    for i,read in enumerate(sys.stdin):
        if i%3_000_000==0:
            print("Reading BAM:", datetime.datetime.now(), i, flush=True)
            
        reads.append(read)

        if len(reads)>=read_batch_n:
            for read in reads:
                parse_sam_line(read)
            reads = []
        
    else:
        for read in reads:
            parse_sam_line(read)
        db.commit()
        del reads

    query_index = None

    # Add database indexes

    db.execute('CREATE INDEX query_name_idx ON query (name);')
    db.execute('CREATE INDEX reference_name_idx ON reference (name);')
    db.execute('CREATE INDEX cluster_name_idx ON cluster (name);')
    db.execute('CREATE INDEX species_name_idx ON species (name);')
    db.execute('CREATE INDEX genome_name_idx ON genome (name);')

    db.execute('CREATE INDEX cds_species_idx ON cds_read_mappings (species);')
    db.execute('CREATE INDEX cds_genome_idx ON cds_read_mappings (genome);')
    db.execute('CREATE INDEX cds_query_idx ON cds_read_mappings (query);')
    db.execute('CREATE INDEX cds_cluster_idx ON cds_read_mappings (cluster);')
    # db.execute('CREATE INDEX cds_reference_idx ON cds_read_mappings (reference);')

    db.execute('CREATE INDEX species_idx ON species_genome_read_mappings (species);')
    db.execute('CREATE INDEX genome_idx ON species_genome_read_mappings (genome);')
    db.execute('CREATE INDEX query_idx ON species_genome_read_mappings (query);')
    # db.execute('CREATE INDEX reference_idx ON species_genome_read_mappings (reference);')


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
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome_list[genome]]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,align_score,ani FROM species_genome_read_mappings WHERE genome={genome};')
        mappings = defaultdict(lambda :(None,None,None,0,0))
        for q,r,rs,re,a,ani in cur:
            if a>=mappings[q][-2]:
                mappings[q] = (r,rs,re,a,ani)
        
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
            cur.execute(f'SELECT query,reference,rstart,rend,align_score,ani FROM species_genome_read_mappings WHERE genome={genome};')
            mappings = defaultdict(set)
            for q,r,rs,re,a,ani in cur:
                mappings[q].add((r,rs,re,a,ani))

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
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome_list[genome]]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,align_score,ani FROM species_genome_read_mappings WHERE genome={genome};')
        mappings = defaultdict(set)
        for q,r,rs,re,a,ani in cur:
            mappings[q].add((r,rs,re,a,ani))
        
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
