import argparse
import gzip
import json

parser = argparse.ArgumentParser(description='Count then number of sequences and number of reads in a fastx file.')
parser.add_argument('-i', "--input_fp", type=str,
                    required=True,
                    help="Input fasta/fastq filepath.")
parser.add_argument('-o', "--output_fp", type=str,
                    required=True,
                    help="Output path for JSON summary stats file.")
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
                yield base_count
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
        yield base_count

if __name__ == '__main__':
    gz = args.input_fp[-3:]=='.gz'
    in_file = gzip.open(args.input_fp, 'rb') if gz else open(args.input_fp, 'rt')
    base_count = 0
    seq_count = 0

    for base_count_ in get_reads(in_file, gz):
        base_count += base_count_
        seq_count += 1

    with open(args.output_fp, 'wt') as f:
        json.dump({'base_count': base_count, 'seq_count': seq_count}, f)


