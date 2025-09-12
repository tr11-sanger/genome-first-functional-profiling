import pysam
import datetime
import argparse

parser = argparse.ArgumentParser(description='Chunk fasta and fastq files by base and reads. If both base and read targets are set then either of the criteria triggers the splitting of a chunk.')
parser.add_argument('-i', "--input_fp", type=str,
                    required=True,
                    help="Input BAM filepath.")
parser.add_argument('-o', "--output_fp", type=str,
                    required=True,
                    help="Output CSV filepath.")
args = parser.parse_args()

def extract_from_bam(bam):
    samfile = pysam.AlignmentFile(bam, "rb")
    for read in samfile:
        if not read.aligned_pairs:
            continue
        query = read.query_name
        reference = read.reference_name
        qstart, rstart = read.aligned_pairs[0]
        qend, rend = read.aligned_pairs[-1]
        ani = 1-(read.get_tag("NM")/read.query_alignment_length)
        alength = qend-qstart
        mapq = read.mapping_quality
        align_score = read.get_tag("AS")
        yield query, reference, rstart, rend, alength, mapq, align_score, ani

if __name__ == '__main__':
    with open(args.output_fp, 'wt') as f:
        f.write("query,reference,rstart,rend,alength,mapq,as\n")
        for i,(q,r,s,e,l,m,s,a) in enumerate(extract_from_bam(args.input_fp)):
            if i%1_000_000==0:
                print(datetime.datetime.now(), i)
            f.write(f"{q},{r},{s},{e},{l},{m},{s},{a}\n")