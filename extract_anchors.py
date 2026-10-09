#!/usr/bin/env python3
import argparse
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract single copy regions (Anchors) from a FASTA file using a canonical k-mer count table."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "-k", "--kmers", required=True, help="Input k-mer table TSV (kmer\\tcount)"
    )
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument(
        "--min-len",
        type=int,
        default=100,
        help="Minimum length for single-copy anchors [default: 100]",
    )
    parser.add_argument(
        "--max-gap",
        type=int,
        default=50,
        help="Maximum gap distance to bridge [default: 50]",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Minimum k-mer count cutoff [default: 1]",
    )
    return parser.parse_args()


def reverse_complement(seq):
    trans = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(trans)[::-1]


def canonical_kmer(seq):
    rc = reverse_complement(seq)
    return seq if seq < rc else rc


def load_kmers(kmers_tsv):
    sys.stderr.write(f"Loading k-mer counts from {kmers_tsv}...\n")
    kmer_counts = {}
    kmer_len = None

    with open(kmers_tsv, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                parts = line.split()
            kmer, count = parts[0], int(parts[1])

            if kmer_len is None:
                kmer_len = len(kmer)

            kmer_counts[kmer] = count

    sys.stderr.write(
        f"Loaded {len(kmer_counts)} k-mers of length {kmer_len}.\n"
    )
    return kmer_counts, kmer_len


def read_fasta(fasta_file):
    current_header = None
    current_seq = []

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    yield current_header, "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            yield current_header, "".join(current_seq)


def process_sequence(
    chrom, seq, kmer_counts, kmer_len, min_count, min_len, max_gap, out_handle
):
    anchor_start = -1
    anchor_end = -1
    in_anchor = False
    num_haptigs = 0
    prev_start = -1

    seq_len = len(seq)

    def write_anchor_if_valid(start, end, haptig_count):
        length = end - start
        if length >= min_len:
            haptig_count += 1
            anchor_seq = seq[start:end]
            header = f">{chrom}_{haptig_count} {chrom}:{start}-{end}\n"
            out_handle.write(header)

            # Write sequence wrapped at 80 characters
            for i in range(0, len(anchor_seq), 80):
                out_handle.write(anchor_seq[i : i + 80] + "\n")
            return 1
        return 0

    # Slide across sequence
    for start in range(seq_len - kmer_len + 1):
        kmer_seq = seq[start : start + kmer_len].upper()

        # Check for non-ATGC characters
        if any(b not in "ACGT" for b in kmer_seq):
            if in_anchor:
                num_haptigs += write_anchor_if_valid(
                    anchor_start, anchor_end, num_haptigs
                )
                in_anchor = False
            prev_start = start
            continue

        # Check k-mer iteration gap (e.g. if previous iteration skipped over 'N's)
        skipped_due_to_n = prev_start != -1 and start != prev_start + 1
        if skipped_due_to_n and in_anchor:
            num_haptigs += write_anchor_if_valid(
                anchor_start, anchor_end, num_haptigs
            )
            in_anchor = False

        canon = canonical_kmer(kmer_seq)
        count = kmer_counts.get(canon, 0)

        is_qualifying = 0 < count <= min_count

        if is_qualifying:
            if not in_anchor:
                # Start new anchor
                anchor_start = start
                anchor_end = start + kmer_len
                in_anchor = True
            else:
                gap_bases = start - anchor_end

                if gap_bases <= max_gap:
                    # Bridge gap and extend anchor end
                    anchor_end = start + kmer_len
                else:
                    # Gap exceeds max_gap: flush active anchor and start new one
                    num_haptigs += write_anchor_if_valid(
                        anchor_start, anchor_end, num_haptigs
                    )
                    anchor_start = start
                    anchor_end = start + kmer_len

        prev_start = start

    # Flush final anchor at end of sequence
    if in_anchor:
        num_haptigs += write_anchor_if_valid(
            anchor_start, anchor_end, num_haptigs
        )

    sys.stderr.write(
        f"Completed processing sequence {chrom}. Total anchors found: {num_haptigs}\n"
    )


def main():
    args = parse_args()

    kmer_counts, kmer_len = load_kmers(args.kmers)

    with open(args.output, "w") as out_handle:
        for chrom, seq in read_fasta(args.input):
            sys.stderr.write(
                f"Processing {chrom} (Length: {len(seq)} bp)...\n"
            )
            process_sequence(
                chrom,
                seq,
                kmer_counts,
                kmer_len,
                args.min_count,
                args.min_len,
                args.max_gap,
                out_handle,
            )


if __name__ == "__main__":
    main()
