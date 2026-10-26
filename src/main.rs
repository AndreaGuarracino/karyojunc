//! karyojunc: band junctions and band k-mer counts from a KaryoScope annotation BED, in one pass.
mod bands;
use clap::Parser;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

/// Band junctions and band k-mer counts from a KaryoScope annotation BED.
#[derive(Parser)]
#[command(version)]
struct Args {
    /// Preset: sr (short reads, -m 20), lr (long reads, -m 100), asm (assemblies, -m 100)
    #[arg(short = 'x', long, value_parser = ["sr", "lr", "asm"], required_unless_present = "min_seg", conflicts_with = "min_seg")]
    preset: Option<String>,
    /// Minimum bp of a band run, instead of -x
    #[arg(short, long)]
    min_seg: Option<u64>,
    /// KaryoScope annotation BED of reads or contigs, plain or gzip, '-' for stdin
    #[arg(short, long, value_name = "BED")]
    input: String,
    /// Reference bands BED (chrom start end label), the bands of the KaryoScope database
    #[arg(short, long, value_name = "BED")]
    bands: String,
    /// Output prefix
    #[arg(short, long, value_name = "PREFIX")]
    output: String,
    /// Read names, only for paired short reads
    #[arg(short, long, value_name = "TXT")]
    names: Option<String>,
}
const NONE: u16 = u16::MAX - 1; // mate without a band
const MULTI: u16 = u16::MAX; // mate with two or more bands

// fast hash for labels, label pairs and hashed read names
#[derive(Default)]
struct Fx(u64);
impl Hasher for Fx {
    fn finish(&self) -> u64 { self.0 }
    fn write(&mut self, b: &[u8]) { for &x in b { self.write_u64(x as u64) } }
    fn write_u64(&mut self, x: u64) { self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x9E3779B97F4A7C15) }
    fn write_u32(&mut self, x: u32) { self.write_u64(x as u64) }
    fn write_usize(&mut self, x: usize) { self.write_u64(x as u64) }
}
type Map<K, V> = HashMap<K, V, BuildHasherDefault<Fx>>;

fn fnv(s: &[u8]) -> u64 { s.iter().fold(0xcbf29ce484222325, |h, &b| (h ^ b as u64).wrapping_mul(0x100000001b3)) }
// unsigned integer as Rust's u64 parse reads it (optional '+')
fn num(x: &[u8]) -> Option<u64> {
    let d = x.strip_prefix(b"+").unwrap_or(x);
    if d.is_empty() { return None; }
    d.iter().try_fold(0u64, |n, &c| if c.is_ascii_digit() { n.checked_mul(10)?.checked_add((c - b'0') as u64) } else { None })
}

// plain text, gzip through bgzip, or '-' for stdin
fn open(p: &str) -> Box<dyn BufRead> {
    if p == "-" { return Box::new(BufReader::with_capacity(1 << 22, std::io::stdin())); }
    let mut f = std::fs::File::open(p).unwrap_or_else(|e| panic!("{p}: {e}"));
    let mut m = [0u8; 2];
    let gz = f.read(&mut m).unwrap_or(0) == 2 && m == [0x1f, 0x8b];
    if !gz { return Box::new(BufReader::with_capacity(1 << 22, std::fs::File::open(p).unwrap())); }
    let mut c = std::process::Command::new("bgzip").args(["-dc", "-@", "4", p]).stdout(std::process::Stdio::piped()).spawn()
        .unwrap_or_else(|e| panic!("bgzip {p}: {e}"));
    let out = c.stdout.take().unwrap();
    Box::new(BufReader::with_capacity(1 << 22, Gz(c, out)))
}

// bgzip output; a failed bgzip is an error at end of stream
struct Gz(std::process::Child, std::process::ChildStdout);
impl std::io::Read for Gz {
    fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
        let n = self.1.read(b)?;
        if n == 0 && !b.is_empty() {
            let s = self.0.wait()?;
            if !s.success() { return Err(std::io::Error::other(format!("bgzip {s}"))); }
        }
        Ok(n)
    }
}

#[derive(Default)]
struct St {
    ids: Map<Vec<u8>, u16>,
    names: Vec<String>,
    bp: Vec<u64>,                     // k-mer positions per label
    loc: Vec<Option<(u8, u64)>>,      // chromosome id and midpoint per label
    chroms: Vec<String>,
    pairs: Map<u32, [u64; 6]>,        // discordant, split, revisit, dom_a, dom_b, inverted
    open: Map<u64, u16>,              // mates waiting for their pair
    reads: u64,
    frags: u64,
}

impl St {
    fn id(&mut self, l: &[u8], b: &bands::Bands) -> u16 {
        if let Some(&i) = self.ids.get(l) { return i; }
        let i = self.names.len() as u16;
        assert!(i < NONE, "more than {NONE} labels");
        let s = String::from_utf8_lossy(l).into_owned();
        let loc = b.at(&s).map(|(c, m)| {
            let c = self.chroms.iter().position(|x| *x == c).unwrap_or_else(|| { self.chroms.push(c); self.chroms.len() - 1 });
            assert!(c < 256, "more than 256 chromosomes"); (c as u8, m)
        });
        self.loc.push(loc);
        self.ids.insert(l.to_vec(), i); self.names.push(s); self.bp.push(0);
        i
    }

    // one read: seq holds its band runs (label, bp)
    fn flush(&mut self, seq: &[(u16, u64)], mate: Option<u64>) {
        let St { names, loc, pairs, open, reads, frags, .. } = self;
        let key = |a: u16, b: u16| if names[a as usize] <= names[b as usize] { (a as u32) << 16 | b as u32 } else { (b as u32) << 16 | a as u32 };
        *reads += 1;
        if seq.len() > 1 {
            let mut adj: Vec<u32> = seq.windows(2).map(|w| key(w[0].0, w[1].0)).collect();
            adj.sort_unstable(); adj.dedup();
            for k in adj { pairs.entry(k).or_default()[1] += 1; }
            // orient on the dominant chromosome (most bp), p to q
            let mut by: Vec<(u8, u64)> = vec![];
            for &(l, n) in seq { if let Some((c, _)) = loc[l as usize] {
                match by.iter_mut().find(|x| x.0 == c) { Some(x) => x.1 += n, None => by.push((c, n)) } } }
            if let Some(&(dom, _)) = by.iter().max_by_key(|x| x.1) {
                let on = |l: u16| loc[l as usize].filter(|x| x.0 == dom).map(|x| x.1 as i64);
                let mut o: Vec<u16> = seq.iter().map(|x| x.0).collect();
                if o.windows(2).filter_map(|w| Some((on(w[1])? - on(w[0])?).signum())).sum::<i64>() < 0 { o.reverse(); }
                for w in o.windows(2) {
                    let k = key(w[0], w[1]);
                    let e = pairs.get_mut(&k).unwrap();
                    e[3] += on((k >> 16) as u16).is_some() as u64;
                    e[4] += on(k as u16).is_some() as u64;
                    e[5] += matches!((on(w[0]), on(w[1])), (Some(x), Some(y)) if y < x) as u64;
                }
            }
            // revisit: back to a band already left (foldback, inversion, tandem duplication)
            for i in 2..seq.len() {
                if seq[..i - 1].iter().any(|x| x.0 == seq[i].0) { pairs.entry(key(seq[i - 1].0, seq[i].0)).or_default()[2] += 1; }
            }
        }
        let Some(h) = mate else { return };
        let me = match seq { [] => NONE, [(l, _)] => *l, _ => MULTI };
        match open.remove(&h) {
            None => { open.insert(h, me); }
            Some(o) => { *frags += 1; if o < NONE && me < NONE { pairs.entry(key(o, me)).or_default()[0] += 1; } }
        }
    }
}

fn main() {
    let args = Args::parse();
    let min_seg = args.min_seg.unwrap_or(if args.preset.as_deref() == Some("sr") { 20 } else { 100 });
    let (out, bands) = (&args.output, bands::Bands::load(&args.bands));
    let mut names = args.names.as_deref().map(open);
    let mut bed = open(&args.input);
    let mut st = St::default();
    let (mut line, mut nm, mut seq) = (vec![], vec![], Vec::<(u16, u64)>::new());
    let (mut cur, mut next_name, mut mate, mut contig, mut nc) = (u64::MAX, 0u64, None, vec![], 0u64);
    loop {
        line.clear();
        let eof = bed.read_until(b'\n', &mut line).expect("read BED") == 0;
        let mut f = line.trim_ascii_end().splitn(5, |&c| c == b'\t');
        let (c0, c1, c2, c3) = (f.next().unwrap(), f.next(), f.next(), f.next());
        assert!(eof || c3.is_some(), "bad BED line: {}", String::from_utf8_lossy(&line));
        // reads are numbered; contigs are counted by name
        let rank = if eof { u64::MAX - 1 } else {
            num(c0).unwrap_or_else(|| { if c0 != contig { contig = c0.to_vec(); nc += 1; } nc })
        };
        if rank != cur {
            if cur != u64::MAX { st.flush(&seq, mate); seq.clear(); }
            if eof { break; }
            assert!(cur == u64::MAX || rank > cur, "ranks not ascending: {rank} after {cur}");
            cur = rank;
            mate = names.as_mut().map(|r| {
                while next_name <= rank {
                    nm.clear();
                    assert!(r.read_until(b'\n', &mut nm).expect("read names") > 0, "names end before rank {rank}");
                    next_name += 1;
                }
                let s = nm.trim_ascii_end();
                fnv(if s.len() > 2 && s[s.len() - 2] == b'/' { &s[..s.len() - 2] } else { s }) // drop /1, /2
            });
        }
        let (c1, c2, c3) = (c1.unwrap(), c2.unwrap(), c3.unwrap());
        let n = num(c2).expect("BED end").wrapping_sub(num(c1).expect("BED start"));
        let i = st.id(c3, &bands);
        st.bp[i as usize] += n;
        if n >= min_seg && bands::is_band(c3) {
            match seq.last_mut() { Some(x) if x.0 == i => x.1 += n, _ => seq.push((i, n)) }
        }
    }
    if st.reads == 0 { eprintln!("karyojunc: no reads in {}", args.input); std::process::exit(2); }

    let tsv = |ext: &str| BufWriter::new(std::fs::File::create(format!("{out}.{ext}.tsv")).unwrap_or_else(|e| panic!("{out}.{ext}.tsv: {e}")));
    let mut w = tsv("band_counts");
    writeln!(w, "label\tkmer_positions").unwrap();
    for (n, b) in st.names.iter().zip(&st.bp) { writeln!(w, "{n}\t{b}").unwrap(); }
    let mut w = tsv("band_pairs");
    writeln!(w, "label_a\tlabel_b\tband_a\tband_b\tedge_class\tsupport\tdiscordant\tsplit\trevisit\tdom_a\tdom_b\tinverted").unwrap();
    let nm = |k: u32| (&st.names[(k >> 16) as usize], &st.names[(k & 0xffff) as usize]);
    let mut v: Vec<_> = st.pairs.iter().collect();
    v.sort_by(|x, y| (y.1[0] + y.1[1]).cmp(&(x.1[0] + x.1[1])).then_with(|| nm(*x.0).cmp(&nm(*y.0))));
    let (mut seen, mut cls) = (HashMap::new(), HashMap::<&str, u64>::new());
    for (&k, c) in v {
        let ((la, lb), s) = (nm(k), c[0] + c[1]);
        let class = bands.classify(la, lb);
        *cls.entry(class).or_default() += s;
        if class == "reference_adjacent" { seen.insert((la.as_str(), lb.as_str()), s); }
        writeln!(w, "{la}\t{lb}\t{}\t{}\t{class}\t{s}\t{}\t{}\t{}\t{}\t{}\t{}", bands::band(la), bands::band(lb), c[0], c[1], c[2], c[3], c[4], c[5]).unwrap();
    }
    let mut w = tsv("adjacency");
    writeln!(w, "label_a\tlabel_b\tchrom\tsupport").unwrap();
    for (la, lb) in &bands.adj {
        let c = la.bytes().take_while(u8::is_ascii_digit).count().clamp(1, 2);
        writeln!(w, "{la}\t{lb}\t{}\t{}", &la[..c], seen.get(&(la.as_str(), lb.as_str())).unwrap_or(&0)).unwrap();
    }
    let mut cl: Vec<_> = cls.into_iter().collect();
    cl.sort_by(|x, y| y.1.cmp(&x.1).then(x.0.cmp(y.0)));
    eprintln!("karyojunc: {} reads, {} mate pairs, {} unpaired, {} labels, {} pairs; support by class: {}",
        st.reads, st.frags, st.open.len(), st.names.len(), st.pairs.len(),
        cl.iter().map(|(k, n)| format!("{k}={n}")).collect::<Vec<_>>().join(" "));
}
