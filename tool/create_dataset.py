#!/usr/bin/env python3

from pathlib import Path
import argparse
import random

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Create a collection from fasta files"
    )

    parser.add_argument('source_dir', help='Directory of the fasta files')
    parser.add_argument("-n", "--nseq", type=int, help="Number of sequence to add")
    parser.add_argument('-d', '--dest_file', help='Location of dest file (default: `data`)', default=None)

    args = parser.parse_args()

    seq_delim = '\30'
    doc_delim = '\0'

    s_dir = Path(args.source_dir)

    # If the destination file wasn't passed, then assume we want to
    # create a new file named 'data'
    d_file = Path(args.dest_file if args.dest_file is not None else 'data')

    with open(d_file, 'w') as writer:
        for filename in s_dir.rglob('*.fasta'):
            print(filename)
            cnt = 0
            with open(filename, 'r') as reader:
                reader.readline()  # Skip heading of first sequence
                for line in reader:
                    if line.startswith('>'):
                        writer.write(seq_delim)
                        cnt += 1
                        if args.nseq and args.nseq <= cnt:
                            break
                        else:
                            continue

                    line = line.rstrip()  # Skip new line characters

                    writer.write(line)

            writer.write(doc_delim)


    # Generate patterns
    end = d_file.stat().st_size - 2
    with open('patterns', 'w') as writer:
        with open(d_file, 'r', errors='ignore') as reader:
            i = 0
            while i < 1000:
                # pattern_size = random.choice([8, 12, 16])
                pattern_size = random.choice([10])
                pos = random.randint(0, end - pattern_size)

                reader.seek(pos)
                pattern = reader.read(pattern_size)

                if pattern.find(seq_delim) == -1 and pattern.find(doc_delim) == -1:
                    writer.write(pattern + '\n')
                    i += 1
