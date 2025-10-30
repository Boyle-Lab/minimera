#!/usr/bin/env python3

import sys
import random
import csv

def avg(xs):
    return sum(xs) / len(xs)

def med(xs):
    xs = sorted(xs)

    if len(xs) % 2 == 0:
        return xs[len(xs)//2]
    else:
        a = len(xs)//2
        b = a+1
        return (xs[a] + xs[b]) / 2

def parse_foldbacks(foldbacks_path):
    with open(foldbacks_path) as f:
        r = csv.reader(f)
        next(r)
        return {id: int(fp) for id, _, kind, _, fp in r if kind == 'foldback'}

def output_foldback_data(read_id, read_qs, fp):
    qscores = [ord(q) - ord('!') for q in read_qs]
    mean_q = avg(qscores)
    median_q = med(qscores)

    for i, q in enumerate(qscores):
        x = i - fp
        print(f"{x},{q},{q - mean_q},{q - median_q}")

def run(fastq_path, foldbacks_path):
    foldbacks = parse_foldbacks(foldbacks_path)
    n=0

    with open(fastq_path) as f:
        while True:
            try:
                read_id = next(f)[1:].split()[0]
            except StopIteration:
                break

            read_seq = next(f).rstrip('\n')
            try:
                next(f)
            except StopIteration:
                print(read_id, file=sys.stderr)
                break
            read_qs = next(f).rstrip('\n')

            if read_id in foldbacks and random.random() < 0.10:
                # print(f'n {n}\n', file=sys.stderr)
                output_foldback_data(read_id, read_qs, foldbacks[read_id])
                n += 1
                if n % 100 == 0:
                    print(n, file=sys.stderr)
                # if n > 2000:
                #     return

_, fastq, foldbacks = sys.argv
run(fastq, foldbacks)
