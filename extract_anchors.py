#!/usr/bin/env python3
import argparse
import sys

# Pre-allocated 256-byte translation table for fast DNA reverse complementation
REV_COMP_TRANS = bytes.maketrans(
    b"ACGTacgtNn",
    b"TGCAtgcaNn"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract single copy regions (Anchors) from a FASTA file using a canonical k-mer count table."
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "-k", "--kmers", required=True, help="Input k-mer table TSV (kmer\tcount)"
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


def reverse_complement(seq_bytes: bytes) -> bytes:
    """Fast reverse complementation operating on raw bytes."""
    return seq_bytes.translate(REV_COMP_TRANS)[::-1]


def canonical_kmer(seq_bytes: bytes) -> bytes:
    rc = reverse_complement(seq_bytes)
    return seq_bytes if seq_bytes < rc else rc


def parse_fasta_bytes(fasta_file):
    """Yields (chrom_bytes, seq_bytearray) reading in binary mode."""
    with open(fasta_file, "rb") as f:
        chrom = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith(b">"):
                if chrom is not None:
                    yield chrom, b"".join(seq_parts)
                chrom = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
        if chrom is not None:
            yield chrom, b"".join(seq_parts)


def main():
    args = parse_args()

    kmers = {}
    k = -1
    with open(args.kmers, "rb") as kmer_file:
        for line in kmer_file:
            parts = line.strip().split(b"\t")
            if len(parts) < 2:
                continue
            kmer, count = parts[0], int(parts[1])
            kmers[kmer] = count
            if k == -1:
                k = len(kmer)
                print(f"Detected k-mer length: {k}", file=sys.stderr)
            elif len(kmer) != k:
                print(f"Error: Inconsistent k-mer lengths in {args.kmers}", file=sys.stderr)
                sys.exit(1)
    print(f"Loaded {len(kmers)} k-mers from {args.kmers}", file=sys.stderr)

    with open(args.output, "wb") as out_handle:
        for chrom, seq in parse_fasta_bytes(args.input):
            print(f"Processing {chrom.decode()}: {len(seq)} bp...", file=sys.stderr)
            anchor = bytearray()
            gap = bytearray()
            anchor_start = 0

            for i in range(len(seq) - k + 1):
                kmer = seq[i : i + k]  # Slice view (bytes)

                if b"N" in kmer:
                    if len(anchor) >= args.min_len:
                        anchor_end = anchor_start + len(anchor)
                        header = f">{chrom.decode()}:{anchor_start}-{anchor_end}\n".encode()
                        out_handle.write(header + anchor + b"\n")
                    anchor.clear()
                    gap.clear()
                    continue

                canon_kmer = canonical_kmer(kmer)
                is_single_copy = (canon_kmer in kmers) and (kmers[canon_kmer] <= args.min_count)

                if is_single_copy:
                    if len(anchor) == 0:
                        # Start a brand new anchor
                        anchor.extend(kmer)
                        anchor_start = i
                        gap.clear()
                    else:
                        if len(gap) > args.max_gap:
                            # Gap too wide; flush existing anchor and start new
                            if len(anchor) >= args.min_len:
                                anchor_end = anchor_start + len(anchor)
                                header = f">{chrom.decode()}:{anchor_start}-{anchor_end}\n".encode()
                                out_handle.write(header + anchor + b"\n")
                            anchor.clear()
                            anchor.extend(kmer)
                            anchor_start = i
                        else:
                            # Bridge gap and append new last byte
                            anchor.extend(gap)
                            anchor.append(kmer[-1])
                        gap.clear()
                else:
                    # Accumulate gap only if currently inside an anchor
                    if len(anchor) > 0:
                        gap.append(kmer[-1])

            # Flush remaining anchor at end of sequence
            if len(anchor) >= args.min_len:
                anchor_end = anchor_start + len(anchor)
                header = f">{chrom.decode()}:{anchor_start}-{anchor_end}\n".encode()
                out_handle.write(header + anchor + b"\n")


if __name__ == "__main__":
    main()