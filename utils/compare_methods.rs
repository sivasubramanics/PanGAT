// Compile with: rustc -C opt-level=3 -C target-cpu=native -C codegen-units=1 compare_methods.rs -o compare_methods

use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::process;
use std::time::Instant;

const HELP: &str = "\
compare_methods - benchmark two BED-based single-copy k-mer methods against k-mer counts

USAGE:
    compare_methods [OPTIONS] <reference.fasta> <kmers_file> <bbmap.bed> <java.bed> [output_file]

ARGUMENTS:
    <reference.fasta>   Reference genome (FASTA, plain text; multi-line records OK)
    <kmers_file>        Text file, one 'KMER COUNT' pair per line (whitespace separated).
                        K-mers are expected to be in canonical form.
    <bbmap.bed>         BED intervals (0-based, half-open) from the BBMap method
    <java.bed>          BED intervals (0-based, half-open) from the Java method
    [output_file]       Per-position TSV table (required unless --no-output is given)

OPTIONS:
    -k, --kmer-size <N>   k-mer length, 1..=32 (default: 31)
    -n, --no-output       Skip the per-position TSV (much faster, far less disk);
                          only the benchmark summary is printed
    -h, --help            Show this help and exit

OUTPUT TSV COLUMNS:
    chromosome  position  kmer  count  is_canonical  bbmap_found  java_found

METRICS (computed on canonical-orientation k-mers only):
    true positive  : count == 1 and position inside a BED interval
    false positive : count != 1 and position inside a BED interval
    recall denominator: all k-mers with count == 1

EXAMPLES:
    compare_methods ref.fa kmers.txt bbmap.bed java.bed out.tsv
    compare_methods -n -k 31 ref.fa kmers.txt bbmap.bed java.bed
";

// ---------- lookup table: ASCII base -> 2-bit code (255 = invalid/N) ----------
const fn make_lut() -> [u8; 256] {
    let mut t = [255u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}
static LUT: [u8; 256] = make_lut();

// ---------- fast integer hasher for packed k-mers ----------
#[derive(Default)]
struct KHasher(u64);

impl Hasher for KHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x517cc1b727220a95);
        }
    }
    #[inline]
    fn write_u64(&mut self, x: u64) {
        let mut z = x.wrapping_add(0x9E3779B97F4A7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        self.0 = z ^ (z >> 31);
    }
}

type KMap = HashMap<u64, u32, BuildHasherDefault<KHasher>>;
type BedMap = HashMap<String, Vec<(u64, u64)>>;

// ---------- helpers ----------
fn die(msg: &str) -> ! {
    eprintln!("error: {}", msg);
    process::exit(1);
}

fn usage_error(msg: &str) -> ! {
    eprintln!("error: {}\n", msg);
    eprintln!("Usage: compare_methods [OPTIONS] <reference.fasta> <kmers_file> <bbmap.bed> <java.bed> [output_file]");
    eprintln!("Try 'compare_methods --help' for more information.");
    process::exit(2);
}

fn open(path: &str) -> File {
    File::open(path).unwrap_or_else(|e| die(&format!("cannot open '{}': {}", path, e)))
}

fn trim(mut s: &[u8]) -> &[u8] {
    while let [f, rest @ ..] = s {
        if f.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    while let [rest @ .., l] = s {
        if l.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    s
}

/// Pack a k-mer string into 2 bits/base (first base = most significant). None if invalid.
fn encode(s: &[u8], k: usize) -> Option<u64> {
    if s.len() != k {
        return None;
    }
    let mut v = 0u64;
    for &b in s {
        let c = LUT[b as usize];
        if c == 255 {
            return None;
        }
        v = (v << 2) | c as u64;
    }
    Some(v)
}

fn push_num(buf: &mut Vec<u8>, mut n: u64) {
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

// ---------- loaders ----------
fn load_kmers(path: &str, k: usize) -> KMap {
    let t0 = Instant::now();
    eprintln!("[kmers] Loading k-mers from {}...", path);
    let f = open(path);
    let est = f.metadata().map(|m| m.len() as usize / (k + 6)).unwrap_or(0);
    let mut map: KMap = HashMap::with_capacity_and_hasher(est, Default::default());
    let mut rdr = BufReader::with_capacity(1 << 20, f);
    let mut line = Vec::new();
    loop {
        line.clear();
        match rdr.read_until(b'\n', &mut line) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => die(&format!("reading '{}': {}", path, e)),
        }
        let mut it = line
            .split(|b| b.is_ascii_whitespace())
            .filter(|s| !s.is_empty());
        let (a, b, c) = (it.next(), it.next(), it.next());
        if c.is_some() {
            continue;
        }
        if let (Some(a), Some(b)) = (a, b) {
            let cnt = std::str::from_utf8(b).ok().and_then(|s| s.parse::<u64>().ok());
            if let (Some(key), Some(cnt)) = (encode(a, k), cnt) {
                map.insert(key, cnt.min(u32::MAX as u64) as u32);
            }
        }
    }
    eprintln!("[kmers] Loaded {} k-mers in {:.2} seconds.", map.len(), t0.elapsed().as_secs_f64());
    map
}

fn load_bed(path: &str) -> BedMap {
    let t0 = Instant::now();
    eprintln!("[bed] Loading BED intervals from {}...", path);
    let mut rdr = BufReader::with_capacity(1 << 20, open(path));
    let mut map: BedMap = HashMap::new();
    let mut line = Vec::new();
    let (mut total, mut skipped) = (0u64, 0u64);
    loop {
        line.clear();
        match rdr.read_until(b'\n', &mut line) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => die(&format!("reading '{}': {}", path, e)),
        }
        let t = trim(&line);
        if t.is_empty() || t[0] == b'#' || t.starts_with(b"track") || t.starts_with(b"browser") {
            continue;
        }
        let mut it = t.split(|b| b.is_ascii_whitespace()).filter(|s| !s.is_empty());
        let chrom = it.next();
        let s = it.next().and_then(|x| std::str::from_utf8(x).ok()).and_then(|x| x.parse::<u64>().ok());
        let e = it.next().and_then(|x| std::str::from_utf8(x).ok()).and_then(|x| x.parse::<u64>().ok());
        match (chrom, s, e) {
            (Some(c), Some(s), Some(e)) if e > s => {
                map.entry(String::from_utf8_lossy(c).into_owned()).or_default().push((s, e));
                total += 1;
            }
            _ => skipped += 1,
        }
    }
    // sort + merge overlapping/adjacent intervals so a single forward pointer suffices
    for v in map.values_mut() {
        v.sort_unstable();
        let mut merged: Vec<(u64, u64)> = Vec::with_capacity(v.len());
        for &(s, e) in v.iter() {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        *v = merged;
    }
    eprintln!(
        "[bed] Loaded {} intervals in {:.2} seconds{}.",
        total,
        t0.elapsed().as_secs_f64(),
        if skipped > 0 { format!(" ({} malformed/empty lines skipped)", skipped) } else { String::new() }
    );
    map
}

// ---------- evaluation ----------
#[derive(Default)]
struct Counters {
    total_single: u64,
    bb_tp: u64,
    bb_fp: u64,
    bb_fn: u64,
    jv_tp: u64,
    jv_fp: u64,
    jv_fn: u64,
}

fn process_record(
    name: &str,
    seq: &[u8],
    k: usize,
    kmers: &KMap,
    bb: &[(u64, u64)],
    jv: &[(u64, u64)],
    mut out: Option<&mut BufWriter<File>>,
    c: &mut Counters,
) -> u64 {
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let shift = (2 * (k - 1)) as u32;
    let (mut fwd, mut rev, mut run) = (0u64, 0u64, 0usize);
    let (mut bi, mut ji) = (0usize, 0usize);
    let mut n_kmers = 0u64;
    let mut line: Vec<u8> = Vec::with_capacity(128);

    for (i, &ch) in seq.iter().enumerate() {
        let code = LUT[ch as usize];
        if code == 255 {
            // N or any non-ACGT base: restart the window
            run = 0;
            fwd = 0;
            rev = 0;
            continue;
        }
        let b = code as u64;
        fwd = ((fwd << 2) | b) & mask;
        rev = (rev >> 2) | ((3 - b) << shift);
        run += 1;
        if run < k {
            continue;
        }
        let pos = (i + 1 - k) as u64;
        n_kmers += 1;

        let is_canon = fwd <= rev;
        let canon = if is_canon { fwd } else { rev };
        let count = kmers.get(&canon).copied().unwrap_or(0);

        while bi < bb.len() && bb[bi].1 <= pos {
            bi += 1;
        }
        while ji < jv.len() && jv[ji].1 <= pos {
            ji += 1;
        }
        let bbf = bi < bb.len() && bb[bi].0 <= pos;
        let jvf = ji < jv.len() && jv[ji].0 <= pos;

        if is_canon {
            if count == 1 {
                c.total_single += 1;
                if bbf { c.bb_tp += 1 } else { c.bb_fn += 1 }
                if jvf { c.jv_tp += 1 } else { c.jv_fn += 1 }
            } else {
                if bbf { c.bb_fp += 1 }
                if jvf { c.jv_fp += 1 }
            }
        }

        if let Some(w) = out.as_deref_mut() {
            line.clear();
            line.extend_from_slice(name.as_bytes());
            line.push(b'\t');
            push_num(&mut line, pos);
            line.push(b'\t');
            line.extend_from_slice(&seq[i + 1 - k..=i]);
            line.push(b'\t');
            push_num(&mut line, count as u64);
            line.push(b'\t');
            line.push(if is_canon { b'1' } else { b'0' });
            line.push(b'\t');
            line.push(if bbf { b'1' } else { b'0' });
            line.push(b'\t');
            line.push(if jvf { b'1' } else { b'0' });
            line.push(b'\n');
            w.write_all(&line).unwrap_or_else(|e| die(&format!("writing output: {}", e)));
        }
    }
    n_kmers
}

fn scores(tp: u64, fp: u64, total_actual: u64) -> (f64, f64, f64) {
    let p = if tp + fp > 0 { tp as f64 / (tp + fp) as f64 * 100.0 } else { 0.0 };
    let r = if total_actual > 0 { tp as f64 / total_actual as f64 * 100.0 } else { 0.0 };
    let f1 = if p + r > 0.0 { 2.0 * p * r / (p + r) } else { 0.0 };
    (p, r, f1)
}

// ---------- CLI ----------
struct Args {
    fasta: String,
    kmers: String,
    bbmap: String,
    java: String,
    output: Option<String>,
    k: usize,
}

fn parse_args() -> Args {
    let argv: Vec<String> = env::args().skip(1).collect();
    let mut pos: Vec<String> = Vec::new();
    let mut k = 31usize;
    let mut no_output = false;
    let mut i = 0;
    while i < argv.len() {
        match argv[i].as_str() {
            "-h" | "--help" => {
                print!("{}", HELP);
                process::exit(0);
            }
            "-k" | "--kmer-size" => {
                i += 1;
                k = argv
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or_else(|| usage_error("--kmer-size needs an integer value"));
            }
            "-n" | "--no-output" => no_output = true,
            s if s.starts_with('-') && s.len() > 1 => usage_error(&format!("unknown option '{}'", s)),
            _ => pos.push(argv[i].clone()),
        }
        i += 1;
    }
    if !(1..=32).contains(&k) {
        usage_error("k-mer size must be between 1 and 32");
    }
    let need = if no_output { 4 } else { 5 };
    if pos.len() < need || pos.len() > 5 {
        usage_error(&format!("expected {} positional arguments, got {}", need, pos.len()));
    }
    Args {
        fasta: pos[0].clone(),
        kmers: pos[1].clone(),
        bbmap: pos[2].clone(),
        java: pos[3].clone(),
        output: if no_output { None } else { pos.get(4).cloned() },
        k,
    }
}

fn main() {
    let args = parse_args();
    let t_all = Instant::now();
    let k = args.k;

    let kmers = load_kmers(&args.kmers, k);
    let bbmap = load_bed(&args.bbmap);
    let java = load_bed(&args.java);
    let empty: Vec<(u64, u64)> = Vec::new();

    let mut out: Option<BufWriter<File>> = args.output.as_ref().map(|p| {
        let f = File::create(p).unwrap_or_else(|e| die(&format!("cannot create '{}': {}", p, e)));
        let mut w = BufWriter::with_capacity(1 << 22, f);
        w.write_all(b"chromosome\tposition\tkmer\tcount\tis_canonical\tbbmap_found\tjava_found\n")
            .unwrap_or_else(|e| die(&format!("writing output: {}", e)));
        w
    });

    eprintln!("[fasta] Starting FASTA scan and evaluation (k={})...", k);
    let mut counters = Counters::default();
    let (mut chrom_count, mut total_kmers) = (0u64, 0u64);

    let mut finish = |name: &str, seq: &[u8]| {
        let t = Instant::now();
        chrom_count += 1;
        let bb = bbmap.get(name).unwrap_or(&empty);
        let jv = java.get(name).unwrap_or(&empty);
        let n = process_record(name, seq, k, &kmers, bb, jv, out.as_mut(), &mut counters);
        total_kmers += n;
        eprintln!("    -> Finished chromosome {} ({} valid k-mers) in {:.2} seconds.", name, n, t.elapsed().as_secs_f64());
    };

    let mut rdr = BufReader::with_capacity(1 << 22, open(&args.fasta));
    let mut line = Vec::new();
    let mut name: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    loop {
        line.clear();
        match rdr.read_until(b'\n', &mut line) {
            Ok(0) => break,
            Ok(_) => {}
            Err(e) => die(&format!("reading '{}': {}", args.fasta, e)),
        }
        let t = trim(&line);
        if t.first() == Some(&b'>') {
            if let Some(nm) = name.take() {
                finish(&nm, &seq);
            }
            let tok = t[1..].split(|b| b.is_ascii_whitespace()).find(|s| !s.is_empty()).unwrap_or(b"");
            name = Some(String::from_utf8_lossy(tok).into_owned());
            seq.clear();
        } else if name.is_some() {
            seq.extend(t.iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(nm) = name.take() {
        finish(&nm, &seq);
    }

    if let Some(w) = out.as_mut() {
        w.flush().unwrap_or_else(|e| die(&format!("flushing output: {}", e)));
    }

    let el = t_all.elapsed().as_secs_f64();
    eprintln!("\n[done] Processed {} total k-mers across {} chromosomes.", total_kmers, chrom_count);
    eprintln!("[done] Total execution time: {:.2} seconds ({:.2} minutes).", el, el / 60.0);

    let (bp, br, bf) = scores(counters.bb_tp, counters.bb_fp, counters.total_single);
    let (jp, jr, jf) = scores(counters.jv_tp, counters.jv_fp, counters.total_single);
    let _ = (counters.bb_fn, counters.jv_fn); // tracked for completeness

    println!("\n================= BENCHMARK RESULTS =================");
    println!("Total true single-copy canonical k-mers: {}", counters.total_single);
    println!("-----------------------------------------------------");
    println!("{:<15} | {:<15} | {:<15}", "Metric", "BBMap Method", "Java Method");
    println!("{}", "-".repeat(51));
    println!("{:<15} | {:<15} | {:.2}% (Are k-mers inside really unique?)", "Precision", format!("{:.2}%", bp), jp);
    println!("{:<15} | {:<15} | {:.2}% (Did we keep all unique k-mers?)", "Recall", format!("{:.2}%", br), jr);
    println!("{:<15} | {:<15} | {:.4}", "F1-Score", format!("{:.4}", bf), jf);
    println!("=====================================================");
}