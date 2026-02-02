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

parser = argparse.ArgumentParser(description='Parse BAM to generate taxonomic and functional profiles. Reassigns multiple mappings with winner-takes-all strategy. Calculates coverage depth and breadth for each species, genome and CDS.')
parser.add_argument('-r', "--refs", type=str,
                    required=True,
                    help="Reference genomes FASTA filepath.")
parser.add_argument('-1', "--reads1", type=str,
                    required=True,
                    help="Raw reads FASTQ filepath.")
parser.add_argument('-2', "--reads2", type=str,
                    default=None,
                    help="Raw reads FASTQ filepath paired ends.")
parser.add_argument('-g', "--genome_species", type=str,
                    required=True,
                    help="Genome to species mapping TSV filepath.")
parser.add_argument('-m', "--genome_contig_mapping", type=str,
                    required=True,
                    help="Genome to contig mapping TSV filepath.")
parser.add_argument('-o', "--out_fp", type=str,
                    required=True,
                    help="Output file path for SQLite database.")
parser.add_argument('-a', "--min_ani", type=float,
                    default=0.95,
                    help="Minimum ANI for read mapping.")
parser.add_argument('-x', "--remove_paired_suffix", type=int,
                    default=1,
                    help="Remove paired suffix of reads (0=No, 1=Yes).")
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

def read_fasta(filehandle, target_seqs=None, split_header=False):
    header = None
    seq_lines = []

    with filehandle as f:
        for l in f:
            l_strip = l.strip()
            if l_strip[0]=='>':
                if not header is None:
                    if (target_seqs is None) or ((target_seqs is not None) and (header in target_seqs)):
                        yield header, ''.join(seq_lines)
                header = l_strip[1:]
                if split_header:
                    header = header.split()[0]
                seq_lines = []
            else:
                if l_strip:
                    seq_lines.append(l_strip)
        else:
            if (target_seqs is None) or ((target_seqs is not None) and (header in target_seqs)):
                yield header, ''.join(seq_lines)

def read_fastq(filehandle, target_seqs=None, split_header=False):
    with filehandle as f:
        for i,l in enumerate(f):
            if i%4==0:
                header = l.strip()[1:]
                if split_header:
                    header = header.split()[0]
            if i%4==1:
                seq = l.strip()
            if i%4==2:
                continue
            if i%4==3:
                q = l.strip()
                if (target_seqs is None) or ((target_seqs is not None) and (header in target_seqs)):
                    yield header, seq, q


if __name__ == '__main__':
    # Set up database
    db_path = args.out_fp  # ':memory:'
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
            name TEXT,
            species INTEGER
        );
    ''')
    db.commit()
    
    db.execute('''
        CREATE TABLE query (
            idx INTEGER PRIMARY KEY,
            name TEXT,
            pair INTEGER,
            length INTEGER
        );
    ''')
    db.commit()
    
    db.execute('''
        CREATE TABLE reference (
            idx INTEGER PRIMARY KEY,
            name TEXT,
            genome INTEGER,
            length INTEGER
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
            ani REAL,
            ani_gapped REAL,
            ani_gapped_fullread REAL
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
            proportion_l REAL,
            ani REAL,
            ani_gapped REAL,
            ani_gapped_fullread REAL
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
            proportion_l REAL,
            ani REAL,
            ani_gapped REAL,
            ani_gapped_fullread REAL
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
            ani REAL,
            ani_gapped REAL,
            ani_gapped_fullread REAL
        );
    ''')
    db.commit()
    
    
    transaction_count = 0
    
    commit_n = 1_000_000
    query_batch_n = 100_000
    
    
    # load data
    
    # Read species and genomes
    genome2species = {}
    with open(args.genome_species, 'rt') as f:
        for l in f:
            k,v = [v.strip() for v in l.split('\t')]
            genome2species[k] = v

    genome_index = {}
    genome_count = 0
    contig2genome = {}
    species_index = {}
    species_count = 0
    
    def add_genome_to_db(genome, db):
        global genome_count
        global genome_index
        global species_count
        global species_index
        global transaction_count
        global commit_n
    
        if genome not in genome2species:
            genome2species[genome] = genome
        species = genome2species[genome]
        if not species in species_index:
            species_index[species] = int(species_count)
            species_count += 1
            db.execute(f'''
                INSERT INTO species (idx, name) 
                VALUES ({species_index[species]}, "{species}");
            ''')
            transaction_count += 1
            if transaction_count%commit_n == 0:
                db.commit()

        genome_index[genome] = int(genome_count)
        genome_count += 1
        db.execute(f'''
            INSERT INTO genome (idx, name, species) 
            VALUES ({genome_index[genome]}, "{genome}", {species_index[species]});
        ''')
        transaction_count += 1
        if transaction_count%commit_n == 0:
            db.commit()
    
    with open(args.genome_contig_mapping, 'rt') as f:
        for l in f:
            a,b = [v.strip() for v in l.split(',')]
            genome = clean_name(a)
            if not genome in genome_index:
                add_genome_to_db(genome, db)
            contig2genome[b] = genome_index[genome]
    genome2species = {genome_index[k]:species_index[v] for k,v in genome2species.items() if k in genome_index}
    
    db.commit()
    transaction_count = 0
    
    
    # Read reference seqs
    reference_index = {}
    reference_count = 0
    
    def add_reference_to_db(reference, length, db):
        global reference_count
        global reference_index
        global transaction_count
        global commit_n
    
        reference_index[reference] = int(reference_count)
        reference_count += 1
        db.execute(f'''
            INSERT INTO reference (idx, name, length, genome) 
            VALUES ({reference_index[reference]}, "{reference}", {length}, {contig2genome[reference]});
        ''')
        transaction_count += 1
        if transaction_count%commit_n == 0:
            db.commit()
            
    with open(args.refs, 'rt') as f:
        for k,v in read_fasta(f, split_header=True):
            add_reference_to_db(k, len(v), db)
    
    db.commit()
    transaction_count = 0
    
    
    # Read FASTQ reads
    query_index = {}
    query_count = 0
    query_length = {}
    
    def add_query_to_db(query, read2, length, db):
        global query_count
        global query_index
        global transaction_count
        global commit_n
    
        query_index[(query,read2)] = int(query_count)
        query_count += 1
        db.execute(f'''
            INSERT INTO query (idx, name, pair, length) 
            VALUES ({query_index[(query,read2)]}, "{query}", {read2}, {length});
        ''')
        transaction_count += 1
        if transaction_count%commit_n == 0:
            db.commit()
    
    read2 = 0
    with gzip.open(args.reads1) as f:
        for k,s,q in read_fastq(f, split_header=True):
            if bool(args.remove_paired_suffix):
                query = k.decode()[:-2]
            else:
                query = k.decode()
            add_query_to_db(query, read2, len(s), db)
            query_length[query_index[(query, read2)]] = len(s)
    if args.reads2 is not None:
        read2 = 1
        with gzip.open(args.reads2) as f:
            for k,s,q in read_fastq(f, split_header=True):
                if bool(args.remove_paired_suffix):
                    query = k.decode()[:-2]
                else:
                    query = k.decode()
                add_query_to_db(query, read2, len(s), db)
                query_length[query_index[(query, read2)]] = len(s)
    
    db.commit()
    transaction_count = 0

    # Read BAM mappings

    nm_suffix = 'nm:i:'
    as_suffix = 'as:i:'
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
        n_matches = 0
        n_exact_matches = 0
        consumes = (False, False)
        v = []
        for c in cigar_str:
            if c in cigar_codes:
                v_int = int(''.join(v))
                if c in {'M','='}:
                    n_matches += v_int
                if c == '=':
                    n_exact_matches += v_int
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
    
        return qalength, ralength, alength, n_matches, n_exact_matches

    def parse_sam_line(read, remove_paired_suffix=False):
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
        global query_batch_n
        global nm_suffix
        global as_suffix
        global contig2genome
        global genome2species

        if read[0]=='@':
            return
        read_ = read.split('\t')
        flags_int = int(read_[1])
        unmapped = bool(flags_int & 4)
        if unmapped:
            return
        if not read_[9]:
            return
        mapq = int(read_[4])
        if mapq==0:
            return
    
        as_tag = None
        nm_tag = None
        for v in read_[10:]:
            if v[:len(nm_suffix)].lower()==nm_suffix:
                nm_tag = int(v[len(nm_suffix):])
                continue
            if v[:len(as_suffix)].lower()==as_suffix:
                as_tag = int(v[len(as_suffix):])
                continue
            
        query = read_[0]
        if remove_paired_suffix:
            paired = True
            read1 = (query[-1]=='1')
            query = query[:-2]
        else:
            paired = bool(flags_int & 1)
            read1 = bool(flags_int & 64)
        read2 = int(not read1) # if paired else None
        cigar = read_[5]
        qa_length, ra_length, alength, n_matches, n_exact_matches = parse_cigar(cigar)
        forward = flags_int & (16 if read1 else 32) 
        reference = read_[2].split()[0]
        rlength = int(read_[7])
        rstart = int(read_[3])
        rend = rstart+ra_length
        ani = 1-(nm_tag/alength)
        ani_gapped = n_matches/alength
        q_length = query_length[query_index[(query,read2)]]
        ani_gapped_fullread = n_matches/max([alength, q_length])
        align_score = as_tag

        if ani_gapped_fullread<args.min_ani:
            return

        genome = contig2genome[reference]
        species = genome2species[genome] if genome in genome2species else 'unknown'

        db.execute(f'''
            INSERT INTO species_genome_read_mappings (
                species, genome, query, reference,
                rstart, rend, ani, ani_gapped, ani_gapped_fullread
            ) 
            VALUES (
                {species}, {genome}, {query_index[(query,read2)]}, {reference_index[reference]}, 
                {rstart}, {rend}, {ani}, {ani_gapped}, {ani_gapped_fullread}
            );
        ''')
        transaction_count += 1
        
        if transaction_count%commit_n == 0:
            db.commit()

    for i,read in enumerate(sys.stdin):
        if i%1_000_000==0:
            print("Reading mappings:", datetime.datetime.now(), i, flush=True)
        parse_sam_line(read, False)  # bool(args.remove_paired_suffix)
    else:
        db.commit()

    n_queries = len(query_index)
    query_index = None


    # Add database indexes

    db.execute('CREATE INDEX query_name_idx ON query (name);')
    db.execute('CREATE INDEX reference_name_idx ON reference (name);')
    db.execute('CREATE INDEX reference_genome_idx ON reference (genome);')
    db.execute('CREATE INDEX cluster_name_idx ON cluster (name);')
    db.execute('CREATE INDEX species_name_idx ON species (name);')
    db.execute('CREATE INDEX genome_name_idx ON genome (name);')
    db.execute('CREATE INDEX genome_species_idx ON genome (species);')

    db.execute('CREATE INDEX species_cds_species_idx ON species_cds_read_mappings (species);')
    db.execute('CREATE INDEX species_cds_genome_idx ON species_cds_read_mappings (genome);')
    db.execute('CREATE INDEX species_cds_query_idx ON species_cds_read_mappings (query);')
    db.execute('CREATE INDEX species_cds_cluster_idx ON species_cds_read_mappings (cluster);')
    db.execute('CREATE INDEX species_cds_ani_idx ON species_cds_read_mappings (ani);')
    db.execute('CREATE INDEX species_cds_ani_gapped_idx ON species_cds_read_mappings (ani_gapped);')
    db.execute('CREATE INDEX species_cds_ani_gapped_fullread_idx ON species_cds_read_mappings (ani_gapped_fullread);')
    # db.execute('CREATE INDEX cds_reference_idx ON cds_read_mappings (reference);')

    db.execute('CREATE INDEX species_idx ON species_genome_read_mappings (species);')
    db.execute('CREATE INDEX genome_idx ON species_genome_read_mappings (genome);')
    db.execute('CREATE INDEX query_idx ON species_genome_read_mappings (query);')
    db.execute('CREATE INDEX ani_idx ON species_genome_read_mappings (ani);')
    db.execute('CREATE INDEX ani_gapped_idx ON species_genome_read_mappings (ani_gapped);')
    db.execute('CREATE INDEX ani_gapped_fullread_idx ON species_genome_read_mappings (ani_gapped_fullread);')
    # db.execute('CREATE INDEX reference_idx ON species_genome_read_mappings (reference);')

    db.execute('CREATE INDEX top_species_idx ON top_species_genome_read_mappings (species);')
    db.execute('CREATE INDEX top_genome_idx ON top_species_genome_read_mappings (genome);')

    db.commit()
    db.close()

    # gzip
    subprocess.call(f"gzip -c {db_path} > {db_path}.gz", shell=True)

    # delete ungzipped
    os.remove(db_path)
