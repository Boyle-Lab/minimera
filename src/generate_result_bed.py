#!/usr/bin/env python3

import sys, csv, pprint

def parse_fastq(in_fastq):
    result = {}
    with open(in_fastq, 'r') as f:
        while True:
            info = f.readline()
            if not info:
                return result

            f.readline()
            f.readline()
            f.readline()

            # >24d1026d-1523-40da-a577-0cdaf2936e06,unknown,chr2:44592758,False
            read_id, read_class, read_loc, read_is_rev = info[1:].rstrip().split(',')
            read_contig, read_pos = read_loc.split(":")
            result[read_id] = (read_contig,
                               int(read_pos),
                               {"True": True, "False": False}[read_is_rev])

def parse_results(result_csv):
    result = {}
    with open(result_csv, 'r') as f:
        r = csv.reader(f)
        next(r) # header
        for row in r:
            id, is_chimera, pos = row
            if is_chimera == "true":
                result[id] = int(pos)
    return result

def compute_foldback_coordinate(foldback_pos, alignment_pos, is_rev):
    if is_rev:
        return alignment_pos - foldback_pos
    else:
        return alignment_pos + foldback_pos

def output_bed(id, foldback_pos, contig, alignment_pos, is_rev):
    pos = compute_foldback_coordinate(foldback_pos, alignment_pos, is_rev)
    print(f"{contig}\t{pos-1}\t{pos+1}\t{id}")


def run(in_fastq, result_csv):
    read_info = parse_fastq(in_fastq)
    chimeras = parse_results(result_csv)

    for id, foldback_pos in chimeras.items():
        contig, alignment_pos, is_rev = read_info[id]
        output_bed(id, foldback_pos, contig, alignment_pos, is_rev)


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2])

