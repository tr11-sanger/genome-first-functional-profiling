import numpy as np
import re
import json
import gzip
from pathlib import Path
from collections import defaultdict
import datetime
import pysam
import copy
import os
import argparse

parser = argparse.ArgumentParser(description='Parse BAM to generate taxonomic and functional profiles. Reassigns multiple mappings with winner-takes-all strategy. Calculates coverage depth and breadth for each species, genome and CDS.')
parser.add_argument('-b', "--bam", type=str,
                    required=True,
                    help="Read mapping BAM filepath.")
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

if __name__ == '__main__':
    # load data

    loc_index_denom = args.loc_index_denom

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
    species2genome = defaultdict(set)
    for k,v in genome2species.items():
        species2genome[v].add(k)
    species2genome = dict(species2genome)

    contig2species = {k:genome2species[clean_name(v)] for k,v in contig2genome.items()}


    contig_lengths = {k:len(v) for k,v in read_fasta(args.refs).items()}


    # Read BAM
    
    samfile = pysam.AlignmentFile(args.bam, "rb")

    read_species_mappings = defaultdict(set)
    species_genome_read_mappings = defaultdict(lambda :defaultdict(lambda :defaultdict(list)))
    query_list = []
    query_index = {}
    reference_list = []
    reference_index = {}
    genome_list = []
    genome_index = {}
    species_list = []
    species_index = {}
    for i,read in enumerate(samfile):
        if i%3_000_000==0:
            print("Reading BAM:", datetime.datetime.now(), i)
            
        if not read.aligned_pairs:
            continue
            
        query = read.query_name
        read2 = (False if read.is_read1 else True) if read.is_paired else None
        forward = read.is_forward
        reference = read.reference_name
        qstart, rstart = read.aligned_pairs[0]
        qend, rend = read.aligned_pairs[-1]
        ani = 1-(read.get_tag("NM")/read.query_alignment_length)
        alength = qend-qstart
        rlength = read.reference_length
        align_score = read.get_tag("AS")
        
        genome = contig2genome[reference]
        species = genome2species[genome]
        
        if not query in query_index:
            query_index[query] = len(query_list)
            query_list.append(query)
        
        if not reference in reference_index:
            reference_index[reference] = len(reference_list)
            reference_list.append(reference)
        
        if not genome in genome_index:
            genome_index[genome] = len(genome_list)
            genome_list.append(genome)
        
        if not species in species_index:
            species_index[species] = len(species_list)
            species_list.append(species)
            
        rlength = rend-rstart
        if all([ani>=args.min_ani, alength<=args.max_align_length, alength>=args.min_align_length]):
            read_species_mappings[(query_index[query],read2)].add(species_index[species])
            species_genome_read_mappings[species_index[species]][genome_index[genome]][(query_index[query],read2)].append((reference_index[reference], rstart, rend, align_score, ani))
    

    # Species profile 

    species_read_counts = {}
    for species,d in species_genome_read_mappings.items():
        species_read_counts[species] = max([len(vs) for k,vs in d.items()])
    
    assigned_species = set()
    print("Assigning reads to species:", datetime.datetime.now(), len(assigned_species), len(species_genome_read_mappings))
    while len(assigned_species)<len(species_genome_read_mappings):
        top_species = max(species_read_counts.items(), key=lambda x:x[1])[0]
        assigned_species.add(top_species)
        del species_read_counts[top_species]
        
        reassigned_reads = {q_ for _,d in species_genome_read_mappings[top_species].items() for (q,_),ts in d.items() for q_ in [(q,True),(q,False),(q,None)]}
        reassigned_species = {s for q in reassigned_reads if q in read_species_mappings for s in read_species_mappings[q]} - assigned_species
        for s in reassigned_species:
            new_dict = {}
            for g,d in species_genome_read_mappings[s].items():
                new_dict[g] = {(q,p):t for (q,p),t in d.items() if q not in reassigned_reads}
            species_genome_read_mappings[s] = dict(new_dict)
            species_read_counts[s] = max([len(vs) for k,vs in new_dict.items()])
            
        print("Assigning reads to species:", datetime.datetime.now(), len(assigned_species), len(species_genome_read_mappings), species_list[top_species])
    
    species_genomes_coverage = defaultdict(dict)
    for species, genome_mappings in species_genome_read_mappings.items():
        print("Calculating species coverages:", datetime.datetime.now(), species_list[species])
        
        for genome, mappings in genome_mappings.items():    
            contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome_list[genome]]}
    
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
            
            species_genomes_coverage[species][genome] = (float(genome_coverage_depth), float(genome_coverage_breadth), float(genome_expected_breadth), float(genome_coverage_breadth/genome_expected_breadth), len(mappings), len({k for k,_ in mappings}))
    
    species_top_genome_coverage = {k:max(list(d.items()), key=lambda x:x[1][1]) for k,d in species_genomes_coverage.items()}


    # Genome profile

    genome_greedy_mappings = defaultdict(dict)
    for species, genome_read_mappings_ in species_genome_read_mappings.items():
        
        genome_read_mappings = copy.deepcopy(genome_read_mappings_)
        
        assigned_genomes = set()
        read_genome_mappings = defaultdict(dict)
        for k,d in genome_read_mappings.items():
            for r,vs in d.items():
                read_genome_mappings[r][k] = vs
        
        genome_read_counts = {k:len(vs) for k,vs in genome_read_mappings.items()}
        while (len(assigned_genomes)<len(genome_read_mappings)):
            top_genome = max(genome_read_counts.items(), key=lambda x:x[1])[0]
            assigned_genomes.add(top_genome)
            del genome_read_counts[top_genome]
            
            print("Assigning reads to genomes:", datetime.datetime.now(), species_list[species], len(assigned_genomes), len(genome_read_mappings), genome_list[top_genome])
    
            reassign_reads = {k_ for k,_ in genome_read_mappings[top_genome] for k_ in [(k,True),(k,False),(k,None)]}
            reassign_genomes = {v for k in reassign_reads for v in read_genome_mappings[k]} - assigned_genomes
            for k in reassign_genomes:
                new_dict = {read:vs for read,vs in genome_read_mappings[k].items() if read not in reassign_reads}
                genome_read_mappings[k] = dict(new_dict)
                genome_read_counts[k] = len(new_dict)
        
        for k,v in genome_read_mappings.items():
            genome_greedy_mappings[k] = v

    
    genomes_coverage = {}
    for genome, mappings in genome_greedy_mappings.items():
        print("Calculating genome coverages:", datetime.datetime.now(), genome2species[genome_list[genome]], genome_list[genome])
        contig_coverage_depth = {v: np.zeros(contig_lengths[v]) for v in genome2contigs[genome_list[genome]]}
    
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
    
        genomes_coverage[genome] = (float(genome_coverage_depth), float(genome_coverage_breadth), float(genome_expected_breadth), float(genome_coverage_breadth/genome_expected_breadth), len(mappings), len({k for k,_ in mappings}))
    

    # CDS profile

    cluster_list = []
    cluster_index = {}
    species_cds_mappings = {}
    for species, genome_mappings in species_genome_read_mappings.items():
        read_cds_mappings = defaultdict(set)
        cds_read_mappings = defaultdict(lambda :defaultdict(set))
        for genome, mappings in genome_mappings.items():    
            # load and index genome CDSs
            if not genome_list[genome] in genome2cds_fp:
                continue
                
            cds_file = []
            with open(genome2cds_fp[genome_list[genome]], 'rt') as f:
                for l in f:
                    name, cluster, contig, _, start, end, direction, data = [v.strip() for v in l.split('\t')]
                    
                    if cluster not in cluster_index:
                        cluster_index[cluster] = len(cluster_list)
                        cluster_list.append(cluster)
                    
                    if contig not in reference_index:
                        reference_index[contig] = len(reference_list)
                        reference_list.append(contig)
                    
                    data_dict = {}
                    for s in data.split(';'):
                        k,v = s.split('=')
                        data_dict[k] = v
                    partial = data_dict['partial']!='00'
                    cds_file.append((cluster_index[cluster], reference_index[contig], int(start), int(end), partial))
                    
            cds_file_index = defaultdict(lambda :defaultdict(set))
            for i,(cluster, contig, start, end, partial) in enumerate(cds_file):
                for j in range(start//loc_index_denom, (end+1)//loc_index_denom):
                    cds_file_index[contig][j].add(i)
    
            # translate mappings to CDSs
            for k,ts in mappings.items():
                for t in ts:
                    if t[0] not in cds_file_index:
                        continue
                    
                    cdss = {v for i in range(t[1]//loc_index_denom,(t[2]+1)//loc_index_denom) for v in cds_file_index[t[0]][i]}
                    
                    for cds_i in cdss:
                        cluster, contig, start, end, partial = cds_file[cds_i]
                        start_ = max([t[1], start])
                        end_ = min([t[2], end])
                        l = end_-start_
                        cds_l = end-start
                        if l>0:
                            read_cds_mappings[k].add((cluster,partial))
                            cds_read_mappings[(cluster,partial)][k].add((t[1],t[2],start_,end_,start,end,l/cds_l,genome,contig))
                        
        print("Assigning reads to CDSs:", datetime.datetime.now(), species_list[species], len(cds_read_mappings))

        # greedy read reassignment
        cds_coverages = {k: sum([max(vs, key=lambda x:x[-2])[-1] for _,vs in d.items()]) for k,d in cds_read_mappings.items()}
        
        assigned_cdss = set()
        while (len(assigned_cdss)<len(cds_read_mappings)):
            top_cds = max(cds_coverages.items(), key=lambda x:(-x[0][1],x[1]))[0]
            assigned_cdss.add(top_cds)
            del cds_coverages[top_cds]
            
            reassign_reads = cds_read_mappings[top_cds]
            reassign_cdss = {v for k in reassign_reads for v in read_cds_mappings[k]} - assigned_cdss
            for k in reassign_cdss:
                new_dict = {}
                for read,vs in cds_read_mappings[k].items():
                    if read in reassign_reads:
                        starts, ends = zip(*[(max([0,cs-rs]),min([re-rs,ce-rs])) for rs,re,_,_,cs,ce,_,_,_ in reassign_reads[read]])
                        min_start = min(starts)
                        max_end = max(ends)
                        vs_ = set()
                        for v in vs:
                            rs,re,_,_,cs,ce,_,_,_ = v
                            start = max([0,cs-rs])
                            end = min([re-rs,ce-rs])
                            if (start>max_end) or (end<min_start):
                                vs_.add(v)
                        if vs_:
                            new_dict[read] = vs_
                    else:        
                        new_dict[read] = vs
                    
                cds_read_mappings[k] = dict(new_dict)
                cds_coverages[k] = sum([max(vs, key=lambda x:x[-2])[-1] for _,vs in cds_read_mappings[k].items()])
        
        species_cds_mappings[species] = cds_read_mappings
    
    species_cds_cluster_coverage = defaultdict(dict)
    for species, cds_read_mappings in species_cds_mappings.items():
        print("Calculating CDS coverages:", datetime.datetime.now(), species_list[species])
        
        for cds_cluster, read_mappings in cds_read_mappings.items():
            cds_coverage = {}
            for _, mappings in read_mappings.items():
                cds_mappings = defaultdict(list)
                for read_start, read_end, align_start, align_end, cds_start, cds_end, l, genome, contig in mappings:
                    cds_mappings[(genome, contig, cds_start)].append((read_start, read_end, align_start, align_end, cds_start, cds_end, l, genome, contig))
                for k,vs in cds_mappings.items():
                    read_start, read_end, align_start, align_end, cds_start, cds_end, l, genome, contig = max(vs, key=lambda x:x[6])
                    if not k in cds_coverage:
                        cds_coverage[k] = np.zeros(cds_end-cds_start+1, dtype=int)
                    cds_coverage[k][align_start-cds_start:align_end-cds_start+1] += 1
            
            cds_coverage_summary = {}
            for k,vs in cds_coverage.items():
                coverage_depth = vs.mean()
                coverage_breadth = (vs>0).mean()
                mean_depth_ = coverage_depth if coverage_depth<700 else 700
                expected_breadth = 1 - (1/(np.log2(1+np.exp(mean_depth_))))
                cds_coverage_summary[k] = (float(coverage_depth), float(coverage_breadth), float(expected_breadth), float(coverage_breadth/expected_breadth), len(vs))
            
            if cds_coverage_summary:
                species_cds_cluster_coverage[species][cds_cluster] = tuple([max(v) for v in zip(*list(cds_coverage_summary.values()))])


    # Mapping statistics
    
    mapping_statistics = {
        'n_mapped_reads': len(read_species_mappings),
        'n_mapped_reads_after_reassign': len({r for _,gs in species_genome_read_mappings.items() for _,rs in gs.items() for r in rs}),
        'n_mapped_read_pairs_after_reassign': len({r for _,gs in species_genome_read_mappings.items() for _,rs in gs.items() for r,_ in rs}),
        'n_mapped_read_pairs': len({k for k,_ in read_species_mappings}),
    }


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
            for (cds,partial),(depth,b,e,r,l) in d.items():
                f.write(f'{species}\t{cluster_list[cds]}\t{1 if partial else 0}\t{depth}\t{b}\t{e}\t{r}\t{l}\n')
    with gzip.open(out_dir / f"{prefix}species_index.txt.gz", 'wt') as f:
        for s in species_list:
            f.write(f'{s}\n')

    with open(out_dir / f"{prefix}mapping_statistics.json", 'wt') as f:
        json.dump(mapping_statistics, f)
