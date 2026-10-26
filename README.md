# karyojunc

karyojunc finds reads that join two chromosome bands, the footprint of translocations, deletions, duplications and inversions.
It reads the band labels that [KaryoScope](https://github.com/barthel-lab/KaryoScope) puts on every read, and it also counts how much of each band was seen, for copy number.
It works on short reads, long reads and assembled contigs.

## Installation

```shell
git clone https://github.com/AndreaGuarracino/karyojunc
cd karyojunc
cargo build --release
```

The program is `target/release/karyojunc`.

## Run

Paired short reads, with the read names file from `karyoscope annotate --query-names-sidecar`, so that the two reads of a pair are joined:

```shell
karyojunc -x sr -i HG002.smoothed.bed.gz -n HG002.query_names.txt.gz -b chm13_cytobands.bed -o HG002
```

Use `-x lr` for long reads and `-x asm` for assemblies.
`-b` is the band BED of the KaryoScope database used for the annotation (`chrom start end label`).

## Output

`HG002.band_pairs.tsv` has one row per pair of bands found joined, the most supported first.
Its main columns:

- `label_a`, `label_b`: the two bands.
- `edge_class`: what the join means. `reference_adjacent`: the bands are neighbours in the reference, a normal join. `same_chrom_distant` and `inter_chrom`: a rearrangement within one chromosome or between two. `par_artifact`, `acro_satellite_artifact`, `subtelomeric_artifact`: joins that similar sequence on different chromosomes can fake.
- `support`: how many reads show the join.

`HG002.band_counts.tsv` gives, for each band, how many read positions carry it: compare it between samples for copy number.

`HG002.adjacency.tsv` lists every pair of neighbouring bands in the reference and how many reads cover it.

## More details

### Options

- `-x sr|lr|asm`: sets `-m` to 20, 100 and 100.
- `-m`: the shortest stretch of one band, in bp, that counts. Use either `-x` or `-m`.
- `-n`: read names, only for paired short reads.

### Evidence for a join

- `split`: one read runs from one band into the other.
- `discordant`: the two reads of a pair land on the two bands.
- `revisit`: a read goes back to a band it already left, as in a foldback, an inversion or a tandem duplication.
- `dom_a`, `dom_b`: reads that lie mostly on the chromosome of `label_a` or of `label_b`; this tells on which chromosome the join sits.
- `inverted`: reads in which the join runs backwards along that chromosome.

Only band labels are joined; other KaryoScope labels are only counted in `band_counts.tsv`.

## Credits

The join classes follow `band_graph.py` by Griffin Azrak.
