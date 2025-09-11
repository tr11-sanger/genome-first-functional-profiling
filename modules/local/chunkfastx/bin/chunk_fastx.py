import argparse
from collections import defaultdict
import gzip
import re
import math

parser = argparse.ArgumentParser(description='Chunk fasta and fastq files by base and reads. If both base and read targets are set then either of the criteria triggers the splitting of a chunk.')
parser.add_argument('-i', "--input_fp", type=str,
                    required=True,
                    help="Input fasta/fastq filepath.")
parser.add_argument('-o', "--output_pattern", type=str,
                    required=True,
                    help="Output pattern (will take directory, basename and extension and add chunk number and paired-end number).")
parser.add_argument('-z', "--output_gzip", action='store_true',
                    help="Gzip the output chunks. If not specified then will only gzip when inputs are gzipped.")
parser.add_argument('-c', "--no_output_gzip", action='store_false', dest='output_gzip',
                    help="Don't gzip the output chunks. If not specified then will only gzip when inputs are gzipped.")
parser.add_argument('-b', "--base_count_target", type=str,
                    default='', required=False,
                    help="Target number of bases in each chunk (rounded down). K, M, G and T suffixes available (e.g. 1.5K)")
parser.add_argument('-r', "--read_count_target", type=str,
                    default='', required=False,
                    help="Target number of reads in each chunk (rounded down). K, M, G and T suffixes available (e.g. 1.5K).")
parser.add_argument('-n', "--read_batch_n", type=int,
                    default=200000,
                    help="Number of FASTX entries handled together in a read/write batch.")
args = parser.parse_args()

def get_reads(infile, gz: bool):
    prev_line = None
    lines = []
    header = None
    base_count = 0
    content = False
    for line in infile:
        line_str = line.decode('utf8') if gz else line
        line_str_strip = line_str.strip()
        if (not prev_line=='+') and (line_str[0] in {'>','@'}):
            if not header is None:
                yield (base_count, '\n'.join(lines))
            lines = []
            header = line_str[1:].split()[0]
            content = False
            base_count = 0
        if line_str_strip == '+':
            content = False
        if ((prev_line is not None) and (prev_line[0] in {'>','@'})) and (not line_str[0] in {'>','@','+'}):
            content = True
        if content:
            base_count += len(line_str_strip)
        lines.append(line_str_strip)
        prev_line = line_str_strip
    else:
        yield (base_count, '\n'.join(lines))

def convert_size(size_str, size_suffix=None):
    if size_suffix is None:
        size_suffix = {v:i for i,v in enumerate(["", "K", "M", "G", "T"])}
    v, s = re.findall(r"^([0-9.,]+)([a-zA-Z]*)?$", size_str)[0]
    s = s[0] if s else ""
    return int(float(re.sub(r"[^0-9.]", '', v)) * math.pow(1000, size_suffix[s]))

if __name__ == '__main__':
    gz = args.input_fp[-3:]=='.gz'
    in_reads = get_reads(gzip.open(args.input_fp, 'rb') if gz else open(args.input_fp, 'rt'), gz)
    base_count = 0
    read_count = 0
    stop_flag = False

    basename, extension = re.findall(r"^(.*?)\.([^.]+(\.gz)?)$", args.output_pattern)[0][:2]
    chunk_n = 1
    out_gz = (args.output_gzip) if (args.output_gzip is not None) else (extension[-3:]=='.gz')
    out_file = None

    size_suffix = {v:i for i,v in enumerate(["", "K", "M", "G", "T"])}
    read_count_target = convert_size(args.read_count_target, size_suffix) if args.read_count_target else False
    base_count_target = convert_size(args.base_count_target, size_suffix) if args.base_count_target else False

    while True:
        read_stack = []
        for _ in range(args.read_batch_n):
            if stop_flag:
                continue
            try:
                base_count_, v = next(in_reads)
            except StopIteration:
                stop_flag = True
                continue
            read_stack.append((base_count_, v))

        # output reads that have been read in all input files
        for base_count_, v in read_stack:
            if out_file is None:
                out_fp = f"{basename}.chunk-{chunk_n}.{extension}"
                out_file = gzip.open(out_fp, 'wb') if out_gz else open(out_fp, 'wt')

            out_str = v+'\n'
            out_file.write(out_str.encode() if out_gz else out_str)
            read_count += 1
            base_count += base_count_
            # check chunking, reset counts and update out_files
            chunk = False
            if base_count_target:
                if base_count>=base_count_target:
                    chunk = True
            if read_count_target:
                if read_count>=read_count_target:
                    chunk = True
            if chunk:
                base_count = 0
                read_count = 0
                chunk_n += 1
                out_file.close()
                out_file = None

        # exit condition
        if stop_flag:
            break

    if not out_file is None:
        out_file.close()

