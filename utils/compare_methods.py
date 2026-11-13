import sys
import time
from collections import defaultdict

def load_kmers(kmers_file):
    """Load canonical k-mers and their counts into a dictionary with byte keys."""
    start_time = time.time()
    kmer_dict = {}
    with open(kmers_file, 'r') as f:
        print(f"[kmers] Loading k-mers from {kmers_file}...")
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                kmer_str, count_str = parts[0], parts[1]
                try:
                    kmer_dict[kmer_str.encode('ascii')] = int(count_str)
                except ValueError:
                    continue
    elapsed = time.time() - start_time
    print(f"[kmers] Loaded {len(kmer_dict)} k-mers in {elapsed:.2f} seconds.")
    return kmer_dict

def load_bed(bed_file):
    """Load BED file into a dictionary with chromosome bytes keys and integer bins."""
    start_time = time.time()
    bed_dict = defaultdict(lambda: defaultdict(list))
    bin_size = 10000
    total_intervals = 0

    with open(bed_file, 'r') as f:
        print(f"[bed] Loading BED intervals from {bed_file}...")
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            chrom = parts[0].encode('ascii')
            start, end = int(parts[1]), int(parts[2])

            start_bin = start // bin_size
            end_bin = (end - 1) // bin_size
            for b in range(start_bin, end_bin + 1):
                bed_dict[chrom][b].append((start, end))
            total_intervals += 1

    elapsed = time.time() - start_time
    print(f"[bed] Loaded {total_intervals} intervals in {elapsed:.2f} seconds.")
    return bed_dict

def parse_fasta(fasta_file):
    """Parse fasta file and yield (id_bytes, sequence_bytes) tuples."""
    with open(fasta_file, 'r') as f:
        seq_id = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    yield (seq_id, b''.join(seq_lines))
                seq_id = line[1:].split()[0].encode('ascii')
                seq_lines = []
            else:
                seq_lines.append(line.upper().encode('ascii'))
        if seq_id is not None:
            yield (seq_id, b''.join(seq_lines))

TRANSTAB = bytes.maketrans(b'ACGTNacgtn', b'TGCANtgcan')

def canonical(kmer_bytes):
    """Return the canonical form of a byte-string k-mer using fast C translation."""
    rev_comp = kmer_bytes.translate(TRANSTAB)[::-1]
    return kmer_bytes if kmer_bytes < rev_comp else rev_comp

def walk_kmers(sequence_bytes, k):
    """Yield (pos, raw_kmer_bytes, canonical_kmer_bytes, is_canonical) tuples."""
    for i in range(len(sequence_bytes) - k + 1):
        kmer = sequence_bytes[i:i+k]
        if b'N' in kmer or b'n' in kmer:
            continue
        canon = canonical(kmer)
        yield (i, kmer, canon, 1 if canon == kmer else 0)

def main():
    if len(sys.argv) != 6:
        print("Usage: python compare_methods.py <reference.fasta> <kmers_file> <bbmap2.bed> <java.bed> <output_file>")
        sys.exit(1)

    overall_start = time.time()

    kmers = load_kmers(sys.argv[2])
    bbmap_bed = load_bed(sys.argv[3])
    java_bed = load_bed(sys.argv[4])
    bin_size = 10000

    buffer_size = 50000
    buffer = []

    # Evaluation metric counters
    bb_tp, bb_fp, bb_fn = 0, 0, 0
    jv_tp, jv_fp, jv_fn = 0, 0, 0
    total_canonical_single = 0

    print(f"[fasta] Starting FASTA scan and evaluation...")
    chrom_count = 0
    total_kmers_processed = 0

    with open(sys.argv[5], 'w', buffering=8192) as out_f:
        out_f.write("chromosome\tposition\tkmer\tcount\tis_canonical\tbbmap_found\tjava_found\n")

        for chrom_bytes, seq_bytes in parse_fasta(sys.argv[1]):
            chrom_start = time.time()
            chrom_count += 1
            seq_kmers = 0

            chrom_str = chrom_bytes.decode('ascii')

            for pos, raw_kmer_bytes, canon_kmer_bytes, is_canonical in walk_kmers(seq_bytes, 31):
                seq_kmers += 1
                kmer_count = kmers.get(canon_kmer_bytes, 0)

                idx = pos // bin_size
                bbmap_intervals = bbmap_bed.get(chrom_bytes, {}).get(idx, [])
                java_intervals = java_bed.get(chrom_bytes, {}).get(idx, [])

                bbmap_found = 1 if any(start <= pos < end for start, end in bbmap_intervals) else 0
                java_found = 1 if any(start <= pos < end for start, end in java_intervals) else 0

                # Accumulate evaluation metrics on-the-fly for canonical k-mers
                if is_canonical == 1:
                    is_true_single = (kmer_count == 1)
                    if is_true_single:
                        total_canonical_single += 1
                        if bbmap_found: bb_tp += 1
                        else: bb_fn += 1
                        if java_found: jv_tp += 1
                        else: jv_fn += 1
                    else:
                        if bbmap_found: bb_fp += 1
                        if java_found: jv_fp += 1

                raw_kmer_str = raw_kmer_bytes.decode('ascii')
                buffer.append(f"{chrom_str}\t{pos}\t{raw_kmer_str}\t{kmer_count}\t{is_canonical}\t{bbmap_found}\t{java_found}\n")

                if len(buffer) >= buffer_size:
                    out_f.writelines(buffer)
                    buffer.clear()

            chrom_elapsed = time.time() - chrom_start
            print(f"    -> Finished chromosome {chrom_str} ({seq_kmers:,} valid k-mers) in {chrom_elapsed:.2f} seconds.")
            total_kmers_processed += seq_kmers

        if buffer:
            out_f.writelines(buffer)
            buffer.clear()

    total_elapsed = time.time() - overall_start
    print(f"\n[done] Successfully processed {total_kmers_processed:,} total k-mers across {chrom_count} chromosomes.")
    print(f"[done] Total execution time: {total_elapsed:.2f} seconds ({total_elapsed / 60:.2f} minutes).")

    # Compute Final Performance Metrics
    def compute_scores(tp, fp, fn, total_actual):
        precision = (tp / (tp + fp)) * 100 if (tp + fp) > 0 else 0
        recall = (tp / total_actual) * 100 if total_actual > 0 else 0
        f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
        return precision, recall, f1

    bb_prec, bb_rec, bb_f1 = compute_scores(bb_tp, bb_fp, bb_fn, total_canonical_single)
    jv_prec, jv_rec, jv_f1 = compute_scores(jv_tp, jv_fp, jv_fn, total_canonical_single)

    print("\n================= BENCHMARK RESULTS =================")
    print(f"Total true single-copy canonical k-mers: {total_canonical_single:,}")
    print("-----------------------------------------------------")
    print(f"{'Metric':<15} | {'BBMap Method':<15} | {'Java Method':<15}")
    print("-" * 51)
    print(f"{'Precision':<15} | {bb_prec:.2f}%          | {jv_prec:.2f}% (Are k-mers inside really unique?)")
    print(f"{'Recall':<15} | {bb_rec:.2f}%          | {jv_rec:.2f}% (Did we keep all unique k-mers?)")
    print(f"{'F1-Score':<15} | {bb_f1:.4f}         | {jv_f1:.4f}")
    print("=====================================================")

if __name__ == "__main__":
    main()
