#!/usr/bin/env python3

import random
import sys
import subprocess

import bamnostic as bs

from pprint import pprint

def send_igv_command(cmd):
    print(cmd)
    subprocess.run(['nc', '-q', '2', '127.0.0.1', '60151'], input=cmd, encoding='ascii')

def igv_go(read_name, contig, start, end):
    buffer = max(int((end - start) * 0.20), 50)
    start = max(0, start-buffer)
    end = end + buffer
    send_igv_command(f'clearSelections alignment')
    send_igv_command(f'goto {contig}:{start}-{end}')
    send_igv_command(f'selectByName alignment {read_name}')
    send_igv_command(f'group selected')

RC = str.maketrans('ACGT', 'TGCA')

def reverse_complement(seq):
    return seq.translate(RC)[::-1]

def write_output(output_file, r, classification):
    seq = r.query_sequence
    qualities = ''.join(bs.core.offset_qual(r.query_qualities))
    pos, rev = compute_raw_position(r)
    contig = r.reference_name

    if r.is_reverse:
        seq = reverse_complement(seq)
        qualities = qualities[::-1]

    print(f'>{r.read_name},{classification},{contig}:{pos},{rev}', file=output_file)
    print(f'{r.query_sequence}', file=output_file)
    print(f'+',                  file=output_file)
    print(f'{qualities}',        file=output_file, flush=True)


def compute_igv_position(r):
    contig = r.reference_name
    start = r.reference_start

    reverse = r.is_reverse
    if not reverse:
        end = start + r.query_length
    else:
        end = r.reference_end - r.query_length

    return contig, start, start + r.query_length

def compute_raw_position(r):
    reverse = r.is_reverse
    cigar = bs.utils.parse_cigar(r.cigarstring)

    if reverse:
        pos = r.reference_end

        (ctag, _), clen = cigar[-1]
        if ctag == 'BAM_CSOFT_CLIP':
            pos += clen
    else:
        pos = r.reference_start

        (ctag, _), clen = cigar[0]
        if ctag == 'BAM_CSOFT_CLIP':
            pos -= clen

    return pos, reverse

def read_input():
    while True:
        s = input('class? ')
        if s == 'q':
            return None
        elif s == 'f':
            return 'foldback'
        elif s == 'n':
            return 'normal'
        elif s == 'l':
            return 'ligation-artifact'
        elif s == 'w':
            return 'weird'
        elif s == 's':
            return 'skip'

def should_check(r):
    return (
        not r.is_unmapped
        and not r.is_supplementary
        and not r.is_secondary
        and r.mapping_quality > 30
    )

def read_queue(input_queue):
    # samtools view data/9b58328c3b631816942cbe400a807241935897e6.sorted.bam | f 1,3,4 | shuf > data/9b58328c3b631816942cbe400a807241935897e6.queue
    with open(input_queue, 'r', encoding='ascii') as f:
        q = []
        for line in f:
            name, contig, pos = line.split(' ')
            if contig != '*':
                pos = int(pos)
                q.append((name, contig, pos))
        random.shuffle(q)
        return q

def run(input_bam, input_queue, output_fastq):
    q = read_queue(input_queue)

    with open(output_fastq, 'w') as output_file:
        with bs.AlignmentFile(input_bam, 'rb') as bam:
            for i, target in enumerate(q):
                rs = bam.fetch(target[1], target[2]-1, target[2]+1)
                for r in rs:
                    if r.read_name == target[0]:
                        if should_check(r):
                            print('.', flush=True, end='')
                            # contig, start, end = compute_igv_position(r)
                            # igv_go(r.read_name, contig, start, end)

                            classification = 'unknown'
                            # classification = read_input()
                            if classification is None:
                                return
                            elif classification == 'skip':
                                continue
                            else:
                                write_output(output_file, r, classification)


if __name__ == '__main__':
    run(sys.argv[1], sys.argv[2], sys.argv[3])
