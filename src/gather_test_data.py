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

def write_output(output_fastq, output_csv, r, classification):
    seq = r.query_sequence
    qualities = ''.join(bs.core.offset_qual(r.query_qualities))

    if r.is_reverse:
        seq = reverse_complement(seq)
        qualities = qualities[::-1]

    print(f'@{r.read_name}', file=output_fastq)
    print(f'{seq}',          file=output_fastq)
    print(f'+',              file=output_fastq)
    print(f'{qualities}',    file=output_fastq, flush=True)

    print(f'{r.read_name},{classification}', file=output_csv, flush=True)


def compute_igv_position(r):
    contig = r.reference_name
    start = r.reference_start

    reverse = r.is_reverse
    if not reverse:
        end = start + r.query_length
    else:
        end = r.reference_end - r.query_length

    return contig, start, start + r.query_length

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

def is_clip(cigar, threshold):
    (ctag, _), clen = cigar
    return ctag == 'BAM_CSOFT_CLIP' and clen >= threshold

def read_has_some_soft_clip(r, threshold):
    cigar = bs.utils.parse_cigar(r.cigarstring)

    if len(cigar) > 0:
        return is_clip(cigar[0], threshold) or is_clip(cigar[-1], threshold)
    else:
        return False

# Reads with no clipping only have a 2.5% chance of getting through.
SOFTCLIP_THRESHOLD = 40
UNCLIPPED_PENALTY = 0.025

def should_check(r):
    if (
        not r.is_unmapped
        and not r.is_supplementary
        and not r.is_secondary
        and r.mapping_quality >= 10
    ):
        if read_has_some_soft_clip(r, SOFTCLIP_THRESHOLD):
            return True
        else:
            if random.random() < UNCLIPPED_PENALTY:
                return True
            else:
                print('Skipping unclipped read.')
                return False
    else:
        print('Skipping non-primary or low quality read.')
        return False

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

def run(input_bam_path, input_queue_path, output_csv_path, output_fastq_path):
    q = read_queue(input_queue_path)
    n = 0

    with open(output_csv_path, 'w') as output_csv:
        with open(output_fastq_path, 'w') as output_fastq:
            with bs.AlignmentFile(input_bam_path, 'rb') as bam:
                for i, target in enumerate(q):
                    rs = bam.fetch(target[1], target[2]-1, target[2]+1)
                    for r in rs:
                        print('.', flush=True, end='')
                        if r.read_name == target[0]:
                            if should_check(r):
                                print('')
                                # print('.', flush=True, end='')
                                contig, start, end = compute_igv_position(r)
                                igv_go(r.read_name, contig, start, end)

                                classification = read_input()
                                if classification is None:
                                    return
                                elif classification == 'skip':
                                    continue
                                else:
                                    write_output(output_fastq, output_csv, r, classification)
                                    n += 1
                                    print(f'{n} recorded')


if __name__ == '__main__':
    (_, input_bam_path, input_queue_path, output_csv_path, output_fastq_path) = sys.argv
    run(input_bam_path, input_queue_path, output_csv_path, output_fastq_path)
