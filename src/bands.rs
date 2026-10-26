//! Reference bands: adjacencies, chromosome ends, and the edge class of a label pair.
use std::collections::{HashMap, HashSet};

type Span = (String, u64, u64); // contig, start, end
type Ends<'a> = (u64, &'a str, u64, &'a str);

pub struct Bands { spans: HashMap<String, Vec<Span>>, pub adj: Vec<(String, String)>, adjset: HashSet<(String, String)>, term: HashSet<String> }

// "1p36.33" -> ("1", 'p', "36.33"); arm u marks a chromosome feature of a personal map
fn parse(l: &str) -> Option<(&str, char, &str)> {
    let d = l.bytes().take_while(u8::is_ascii_digit).count();
    let d = match d { 0 if l.starts_with(['X', 'Y']) => 1, 0 | 3.. => return None, d => d };
    let arm = l[d..].chars().next()?;
    matches!(arm, 'p' | 'q' | 'u').then(|| (&l[..d], arm, &l[d + 1..]))
}
pub fn is_band(l: &[u8]) -> bool {
    let d = l.iter().take_while(|c| c.is_ascii_digit()).count();
    let d = if d == 0 && matches!(l.first(), Some(b'X' | b'Y')) { 1 } else { d };
    d > 0 && matches!(l.get(d), Some(b'p' | b'q' | b'u'))
}
// parent band: "1p36.33" -> "1p36"
pub fn band(l: &str) -> String { parse(l).map_or(l.into(), |(c, a, r)| format!("{c}{a}{}", r.split('.').next().unwrap())) }
fn chrom(l: &str) -> Option<&str> {
    let d = l.bytes().take_while(u8::is_ascii_digit).count().min(2);
    if d > 0 { Some(&l[..d]) } else if l.starts_with(['X', 'Y']) { Some(&l[..1]) } else { None }
}
fn key(a: &str, b: &str) -> (String, String) { if a <= b { (a.into(), b.into()) } else { (b.into(), a.into()) } }

impl Bands {
    pub fn load(path: &str) -> Bands {
        let txt = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{path}: {e}"));
        let mut spans = HashMap::<String, Vec<Span>>::new();
        for l in txt.lines().filter(|l| !l.trim().is_empty()) {
            let f: Vec<&str> = l.split('\t').collect();
            assert!(f.len() >= 4, "{path}: fewer than 4 columns: {l}");
            let (c, s, e, n) = (f[0], f[1].parse::<u64>().unwrap(), f[2].parse::<u64>().unwrap(), f[3].trim());
            if parse(n).is_none() { continue; }
            for k in [n.to_string(), band(n)] { // one span per label and contig
                let v = spans.entry(k).or_default();
                match v.iter_mut().find(|x| x.0 == c) { Some(o) => { o.1 = o.1.min(s); o.2 = o.2.max(e); } None => v.push((c.into(), s, e)) }
            }
        }
        let (mut adj, mut term) = (vec![], HashSet::new());
        for parent in [false, true] { // sub-bands, then parent bands
            let mut by: HashMap<&str, Vec<(u64, u64, &str)>> = HashMap::new();
            // first and last band of each arm: (start, label, end, label)
            let mut ends: HashMap<(&str, &str, char), Ends> = HashMap::new();
            for (n, v) in spans.iter().filter(|(n, _)| (**n == band(n)) == parent) {
                let (lc, arm, _) = parse(n).unwrap();
                for (c, s, e) in v {
                    by.entry(c).or_default().push((*s, *e, n));
                    let x = ends.entry((c, lc, arm)).or_insert((u64::MAX, n, 0, n));
                    if *s < x.0 { (x.0, x.1) = (*s, n); }
                    if *e > x.2 { (x.2, x.3) = (*e, n); }
                }
            }
            for it in by.values_mut() { it.sort(); adj.extend(it.windows(2).filter(|w| w[0].2 != w[1].2).map(|w| key(w[0].2, w[1].2))); }
            term.extend(ends.into_iter().map(|((_, _, arm), v)| if arm == 'p' { v.1 } else { v.3 }.to_string()));
        }
        adj.sort(); adj.dedup();
        let adjset = adj.iter().cloned().collect();
        Bands { spans, adj, adjset, term }
    }
    // chromosome and midpoint of a label, to orient reads
    pub fn at(&self, l: &str) -> Option<(String, u64)> {
        self.spans.get(l).and_then(|v| v.first()).map(|(_, s, e)| (chrom(&band(l)).unwrap_or(l).to_string(), (s + e) / 2))
    }
    pub fn classify(&self, a: &str, b: &str) -> &'static str {
        let (ba, bb) = (band(a), band(b));
        let (Some(ca), Some(cb)) = (chrom(&ba), chrom(&bb)) else { return "unparsed" };
        if a == b { return "same_band"; }
        let shared = |f: &dyn Fn(i64) -> bool| ca == cb && self.spans.get(a).zip(self.spans.get(b)).is_some_and(|(va, vb)|
            va.iter().any(|x| vb.iter().any(|y| x.0 == y.0 && f(x.1.max(y.1) as i64 - x.2.min(y.2) as i64))));
        if shared(&|gap| gap < 0) { return "overlapping"; }
        if self.adjset.contains(&key(a, b)) || shared(&|gap| (0..100_000).contains(&gap)) { return "reference_adjacent"; }
        const XY: [&str; 6] = ["Xp22", "Yp11", "Xq28", "Yq12", "Xq21", "Yq11"];
        if matches!((ca, cb), ("X", "Y") | ("Y", "X")) && (XY.contains(&ba.as_str()) || XY.contains(&bb.as_str())) { return "par_artifact"; }
        let acro = |x: &str| parse(x).is_some_and(|(c, arm, r)| ["13", "14", "15", "21", "22"].contains(&c) && (arm == 'p' || r.split('.').next() == Some("11")));
        if acro(&ba) && acro(&bb) { return "acro_satellite_artifact"; }
        if self.term.contains(a) && self.term.contains(b) { return "subtelomeric_artifact"; }
        if ca == cb { "same_chrom_distant" } else { "inter_chrom" }
    }
}
