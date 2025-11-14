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
parser.add_argument('-t', "--min_coverage_ratio", type=float,
                    default=0.95,
                    help="Minimum observed:expected coverage breadth ratio to exclude species and reassign their reads.")
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
            ani REAL
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE species_cds_read_mappings (
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
            ani REAL
        );
    ''')
    db.commit()

    db.execute('''
        CREATE TABLE top_species_genome_read_mappings (
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

    n_queries = len(query_index)
    query_index = None

    # Add database indexes

    db.execute('CREATE INDEX query_name_idx ON query (name);')
    db.execute('CREATE INDEX reference_name_idx ON reference (name);')
    db.execute('CREATE INDEX cluster_name_idx ON cluster (name);')
    db.execute('CREATE INDEX species_name_idx ON species (name);')
    db.execute('CREATE INDEX genome_name_idx ON genome (name);')

    db.execute('CREATE INDEX species_cds_species_idx ON species_cds_read_mappings (species);')
    db.execute('CREATE INDEX species_cds_genome_idx ON species_cds_read_mappings (genome);')
    db.execute('CREATE INDEX species_cds_query_idx ON species_cds_read_mappings (query);')
    db.execute('CREATE INDEX species_cds_cluster_idx ON species_cds_read_mappings (cluster);')
    db.execute('CREATE INDEX species_cds_ani_idx ON species_cds_read_mappings (ani);')
    # db.execute('CREATE INDEX cds_reference_idx ON cds_read_mappings (reference);')

    db.execute('CREATE INDEX species_idx ON species_genome_read_mappings (species);')
    db.execute('CREATE INDEX genome_idx ON species_genome_read_mappings (genome);')
    db.execute('CREATE INDEX query_idx ON species_genome_read_mappings (query);')
    db.execute('CREATE INDEX ani_idx ON species_genome_read_mappings (ani);')
    # db.execute('CREATE INDEX reference_idx ON species_genome_read_mappings (reference);')

    db.execute('CREATE INDEX top_species_idx ON species_genome_read_mappings (species);')
    db.execute('CREATE INDEX top_genome_idx ON species_genome_read_mappings (genome);')

    # Species profile 
    print('Assigning reads to species:', datetime.datetime.now(), flush=True)
    all_assigned = False
    while not all_assigned:
        # get top read match per genome
        print('Trimming mappings to top match per genome:', datetime.datetime.now(), len(assigned_species), flush=True)
        cur = db.cursor()
        cur.execute(f'''
            SELECT idx,query,species
            FROM species_genome_read_mappings
            GROUP BY query, genome
            HAVING ROWID = MIN(ROWID)
            ORDER BY ani DESC;'''
        )

        # identify mappings that map to the top species
        print('Trimming mappings to only match to top species:', datetime.datetime.now(), len(assigned_species), flush=True)
        assigned_species = np.ones(n_queries, dtype=int) * -1
        good_mappings = []
        for i,q,s in cur:
            if assigned_species[q]==-1:
                assigned_species[q] = s
                good_mappings.append(i)
            else:
                if s == assigned_species[q]:
                    good_mappings.append(i)

        # populate top-mappings table
        for good_mappings_  in batch(list(good_mappings), query_batch_n):
            db.execute(f'''
                INSERT INTO top_species_genome_read_mappings
                SELECT * FROM species_genome_read_mappings
                WHERE idx IN ({",".join([str(v) for v in good_mappings_])});
            ''')
            transaction_count += len(good_mappings_)
    
            if transaction_count%commit_n == 0:
                db.commit()
        else:
            db.commit()
    
        # calculate genome coverages
        species_genomes_coverage = defaultdict(dict)
        cur = db.cursor()
        cur.execute('SELECT DISTINCT species,genome FROM top_species_genome_read_mappings;')
        for species, genome in cur:
            print('Calculating species coverage:', datetime.datetime.now(), species_list[species], genome_list[genome], flush=True)
            
            contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome_list[genome]]}
            
            cur = db.cursor()
            cur.execute(f'SELECT query,reference,rstart,rend,align_score,ani FROM top_species_genome_read_mappings WHERE genome={genome};')
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

        # remove bottom species, or all species that fail
        present_species = {k for k,(d,b,e,r,n,n_) in species_top_genome_coverage.items() if r>=args.min_coverage_ratio}
    
        print(f'Removing {len(species_top_genome_coverage) - len(present_species)} species that are below coverage threshold:', datetime.datetime.now(), flush=True)

        if len(present_species) < len(species_top_genome_coverage):
            # delete mappings not in present species
            cur = db.cursor()
            cur.execute(f'DELETE FROM top_species_genome_read_mappings;')
            cur.execute(f'DELETE FROM species_genome_read_mappings WHERE species NOT IN ({",".join(present_species)});')
        else:
            all_assigned = True
            # remove non-top read mapping data
            cur.execute(f'DELETE FROM species_genome_read_mappings;')
        db.commit()

    # genome coverage
    genomes_coverage = {}
    cur = db.cursor()
    cur.execute('SELECT DISTINCT species,genome FROM top_species_genome_read_mappings;')
    species_genome = [(s,g) for s,g in cur]
    for species,genome in species_genome:
        print('Calculating genome coverage:', datetime.datetime.now(), species_list[species], genome_list[genome], flush=True)
        
        contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome_list[genome]]}
        
        cur = db.cursor()
        cur.execute(f'SELECT query,reference,rstart,rend,align_score,ani FROM top_species_genome_read_mappings WHERE genome={genome};')
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

    # CDS profile
    species_cds_cluster_coverage = defaultdict(dict)

    cur = db.cursor()
    cur.execute('SELECT species,genome FROM top_species_genome_read_mappings;')
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
            cur.execute(f'SELECT query,reference,rstart,rend,align_score,ani FROM top_species_genome_read_mappings WHERE genome={genome};')
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
                                INSERT INTO species_cds_read_mappings (
                                    species, cluster, query, partial, genome, reference,
                                    rstart, rend, qstart, qend, cstart, cend, proportion_l, ani
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
        cur.execute(f'SELECT cluster,partial,query,genome,reference,proportion_l,rstart,rend,cstart,cend FROM species_cds_read_mappings ORDER BY query, ani DESC;')
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

            if not (g,r,cs) in cds_coverages:
                cds_coverages[(g,r,cs)] = np.zeros(ce-cs)
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
            b, e, l_ = max(zip(bs, es, ls), key=lambda a,_b,_c:a)
            species_cds_cluster_coverage[species][cluster] = (sum(ds), b, e, b/e, l_)

        # move rows over to cds_read_mappings, delete all from species_cds_read_mappings
        # db.execute(f'''
        #     INSERT INTO cds_read_mappings
        #     SELECT * FROM species_cds_read_mappings; 
        # ''')
        # db.commit()
        db.execute(f'''
            DELETE FROM species_cds_read_mappings; 
        ''')
        db.commit()

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
