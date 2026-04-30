# GenomeHandler

![Python](https://img.shields.io/badge/Python-3.11%2B-blue?logo=python&logoColor=white)
![License](https://img.shields.io/badge/License-Apache%202.0-green)
![Elytron](https://img.shields.io/badge/Elytron-Biotech-4A90D9?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAA4AAAAOCAYAAAAfSC3RAAABVklEQVR42rWSvUsDQRDF3+7cR4xEC1FblQQ1QcEmlgeCGKwUCbZa21iIrYV2VhbaCGIp3D+gheDVVioGQbBKobGwiHh6ud2xiLlcEhAFfbAwD/a383Z2gX+Q+Fq/AIpFily8/pHShR44jhHr3hGntU4XUlaC9gBkALY06+Pw9vQAgASgG5tlhDkOAWDTFtvCSq5o6HVovSre98Yn5+uQ83YTdDzwnpXntK1YFJozjHJI669bwjJi/XDK6ITbAyBRRWsdg3CuVS6F4LWwLhuv2MTdLMMgMF6E6yHQ6ZDJrnDOryqlfInAAQ8T7UPRwBge3RuiKVhBkq82DaWGeohuMmfIXuRRG7Ah+tGIEUxSyWmwcwqC541DExAaVMwkrK/vETCnKFK+KSf7x+BLQl43PI+qbGFvip8yyA5EqrgDpK4SyUSvpAfeDXfUHb9P/mP3/m4OG4+ASGkedKBZx+RAAAAAElFTkSuQmCC)

**Minimal, fast, and clean Python library for genomes and gene annotations.**

`genomehandler` loads genomes from FASTA, parses annotations from GFF3 or GenBank, and answers the two questions that come up all day in strain analysis: *"what gene sits here?"* and *"what's the cDNA for this transcript?"* It is designed around biotech workloads at [Elytron Biotech](https://github.com/ElytronBiotech) — microbial strains, prokaryotic and fungal genomes — but works on any organism.

---

## Features

- **Inputs:** FASTA (`.fa`, `.fasta`, `.fna`, optionally `.gz`), GFF3 (both hierarchical `gene/mRNA/exon/CDS` and gene-less Bakta/Prokka), GenBank (`.gbk`, `.gb`, optionally `.gz`).
- **Three core classes:** `GenomeDB` (queries + I/O), `Genome` (memory-efficient sequence store), `Gene` (structure, exons, introns, sequences).
- **Coordinates:** 1-based, inclusive — matches the GFF3/GenBank convention so there is no mental arithmetic at the boundary.
- **Storage:** 1 byte/base via NumPy (`np.uint8`). A 3 Gb genome fits in ~3 GB of RAM, and regions are returned by zero-copy slicing.
- **Queries:** per-chromosome sorted interval arrays + `np.searchsorted` ⇒ **O(log n)** overlap queries; name lookups via a multi-map.
- **Strand-aware slicing:** `get_region(..., strand="-")` returns the reverse complement.
- **No database required:** everything lives in memory, `pip install` and go.

---

## Installation

### From GitHub

```bash
# latest main
pip install "genomehandler @ git+https://github.com/ElytronBiotech/genomehandler.git@main"

# tagged release
pip install "genomehandler @ git+https://github.com/ElytronBiotech/genomehandler.git@v0.1.0"
```

### Editable install for development

```bash
python -m venv .venv && source .venv/bin/activate
pip install -U pip
pip install -e .
```

### With `uv`

```bash
uv pip install -e .
```

**Requires** Python ≥ 3.11, `numpy ≥ 2.2.6`, `biopython ≥ 1.85`.

---

## Quickstart

```python
from genomehandler import GenomeDB

db = GenomeDB()

# 1) Load a genome (FASTA or FASTA.gz)
db.load_genome("R000044.fasta.gz")

# 2) Load annotations — works on both hierarchical GFF3 and Bakta/Prokka "gene-less" GFF3
db.load_gff3("R000044.gff3", prefer="cds")   # or prefer="exon"

# 3) Name lookup (multi-mapping -> list)
for g in db.get_genes_by_name("tmk"):
    print(g.region(), g.strand, g.product)

# 4) Position / range queries
hits = db.genes_at("scaffold_1", 5_600)
span = db.genes_overlapping("scaffold_1", 5_200, 10_500)

# 5) cDNA extraction (strand-aware)
if hits:
    cdna = hits[0].extract_cdna_sequence(db.genome)
    print(f"{hits[0].gene_id}: {len(cdna)} bp")
```

---

## Core concepts

### Coordinates
All public methods are **1-based, inclusive**: `start=1, end=3` selects three bases. This matches the conventions used by GFF3, GenBank, BED (reported), and most genome browsers.

### Names vs. IDs
A gene has two handles:
- **`gene_id`** — the stable, unique key (GFF3 `ID`, `locus_tag`, or an auto-generated string). Use `get_genes_by_id()` for O(1) exact lookup.
- **`gene_name`** — the symbol/label that may not be unique (e.g. two paralogs both named `tmk`). Use `get_genes_by_name()` — it always returns a list.

The same name-index is populated from `Name`, `gene`, and `locus_tag` attributes, so common aliases work.

### Interval index
After `load_gff3` / `load_genbank`, `GenomeDB` builds per-chromosome arrays (`starts`, `ends`, `gids`) sorted by start. Overlap queries use `np.searchsorted` to find candidates in O(log n), then a vectorized mask for the end condition.

### Strand handling
The only place strand matters for sequence extraction is `Genome.get_region(chrom, start, end, strand=...)`. Pass `strand="-"` and you get the reverse complement back — no manual complementing at the call site. `Gene.extract_cdna_sequence()` already passes the gene's strand through for you.

---

## API Reference

### `GenomeDB`

Main entry point. Holds one `Genome`, a `dict[gene_id, Gene]`, a name multi-map, protein tables, and per-chromosome interval indices.

| Method | Returns | Purpose |
|---|---|---|
| `load_genome(fasta_file)` | `list[str]` | Load a FASTA (optionally `.gz`). Returns chromosome names. Raises `FileNotFoundError`. |
| `load_gff3(gff_path, *, prefer="cds")` | `int` | Parse a GFF3. Returns number of genes loaded. `prefer` selects `"cds"` or `"exon"` blocks when both exist. Raises `ValueError` for bad `prefer`, `FileNotFoundError` for missing file. |
| `load_genbank(gbk_file, *, load_genome_if_present=True)` | `int` | Parse a GenBank. Attaches protein sequences/IDs to genes, creates genes from CDS features if `self.genes` is empty, and auto-loads the genome from the GBK `ORIGIN` section when `load_genome_if_present=True` and no genome is yet loaded. |
| `get_genes_by_id(gene_id)` | `list[Gene]` | Exact ID match. Returns `[]` if unknown. |
| `get_genes_by_name(name)` | `list[Gene]` | Name/symbol/locus_tag multi-mapping. Returns `[]` if unknown. |
| `genes_at(chrom, pos)` | `list[Gene]` | All genes covering a 1-based position. |
| `genes_overlapping(chrom, start, end)` | `list[Gene]` | All genes intersecting `[start, end]` (swaps if `start > end`). |
| `describe_position(chrom, pos)` | `list[dict]` | JSON-ready summary per overlapping gene (id, name, location, strand, exons, introns, cDNA, protein_id). Populates `cdna_sequence` on the fly when the genome is loaded. |
| `genes_to_matrix(with_strand=False)` | `(np.ndarray, list[str])` | Export all genes as a 2-D object array + column names for DataFrame construction. |

**Attributes** (read-only access is fine):

| Attribute | Type | Contents |
|---|---|---|
| `genome` | `Genome` | The sequence store. |
| `genes` | `dict[str, Gene]` | `gene_id` → `Gene`. |
| `genes_by_name` | `dict[str, list[str]]` | `name` → list of `gene_id`. |
| `proteins` | `dict[str, dict]` | `protein_id` → `{id, name, gene_name, sequence}`. |
| `protein_by_name` | `dict[str, dict]` | `product name` → same dict. |

---

### `Genome`

Memory-efficient per-chromosome sequence store.

```python
from genomehandler import Genome

g = Genome(storage="bytes")   # default; use "int" for integer-coded storage
```

| Method | Returns | Purpose |
|---|---|---|
| `load_fasta(fasta_file)` | `list[str]` | Load FASTA (optionally `.gz`). Returns chromosome names. |
| `load_sequences(seqs)` | `list[str]` | Bulk-load from `list[(name, sequence)]` (used internally when loading from GenBank). |
| `get_base(chrom, pos)` | `str \| None` | Single base at a 1-based position. Returns `None` for unknown chrom or out-of-range pos. |
| `get_region(chrom, start, end, *, strand="+")` | `str \| None` | Substring for `[start, end]`. Returns `None` if chrom unknown; raises `ValueError` if out of bounds. `strand="-"` → reverse complement. |
| `has_chromosome(chrom)` | `bool` | Membership test. |
| `get_chromosome(chrom)` | `np.ndarray \| None` | Raw `uint8` view — zero-copy; useful for vectorized scanning. |
| `get_chromosome_names()` | `list[str]` | All loaded names. |
| `get_chromosome_size(chrom)` | `int \| None` | Length in bases. |
| `get_genome_size()` | `int` | Total length across all chromosomes. |
| `__str__()` | `str` | Human-readable summary with per-chrom sizes. |

**Storage modes:**
- `"bytes"` (default) — ASCII `uint8`. Full-fidelity case preservation and reverse-complement.
- `"int"` — A/C/G/T/N mapped to 1..5. Same 1 byte/base, slightly faster masking workloads, but `get_region(..., strand="-")` only reverses (no complementing). Prefer `"bytes"` unless you specifically need integer-coded arrays.

---

### `Gene`

Represents one gene plus its transcripts, exons, introns, and optional sequences/protein info.

**Fields**

| Field | Type | Notes |
|---|---|---|
| `gene_id` | `str` | Stable key. |
| `gene_name` | `str` | Display symbol. |
| `chromosome` | `str` | Seqid as given in the annotation. |
| `start` / `end` | `int \| None` | 1-based, inclusive. |
| `strand` | `str \| None` | `"+"`, `"-"`, or `"."`. |
| `transcripts` | `list[str]` | Transcript IDs under this gene. |
| `exons` | `dict[(int, int), list[str]]` | Interval → transcripts that include it. |
| `introns` | `dict[(int, int), list[str]]` | Populated by `calculate_introns()`. |
| `cdna_sequence` | `str \| None` | Populated by `extract_cdna_sequence()` when called without a transcript arg. |
| `protein_id` / `protein_sequence` / `product` | `str \| None` | Filled by `load_genbank()` / GFF3 `product=`. |

**Methods**

| Method | Purpose |
|---|---|
| `add_exon(start, end, isoform=None)` | Add an exon block (defaults isoform to `gene_id`). Keeps the transcripts list unique. |
| `calculate_introns()` | Derive introns from exon gaps per transcript. |
| `set_gene_start_end()` | Tighten `start`/`end` to the span of known exons. |
| `extract_cdna_sequence(genome, transcript=None)` | Concatenate exons (strand-aware). With no transcript, returns the full union and stores it in `cdna_sequence`. With a transcript, returns just that isoform and leaves the cached field alone. |
| `as_tuple(with_strand=False)` | `(gene_id, gene_name, chromosome, start, end[, strand])`. |
| `region()` | `"chrom:start-end"` string. |
| `exons_as_tuples(sort=True)` / `introns_as_tuples(sort=True)` | Flat interval lists. |
| `__call__()` | Shortcut for `as_tuple()`. |
| `__str__()` | `"Gene <name> (<id>): chrom:start-end (<strand>)"`. |

---

## Worked examples

### 1. Eukaryotic genome with a hierarchical GFF3

```python
from genomehandler import GenomeDB

db = GenomeDB()
db.load_genome("fungus.fasta.gz")
n = db.load_gff3("fungus.gff3", prefer="exon")   # exon-aware for UTR analysis
print(f"Loaded {n} genes across {len(db.genome.get_chromosome_names())} contigs")

# Per-transcript cDNA
g = db.get_genes_by_name("ACT1")[0]
for tid in g.transcripts:
    seq = g.extract_cdna_sequence(db.genome, transcript=tid)
    print(tid, len(seq))
```

### 2. Bacterial annotation (Bakta / Prokka gene-less GFF3)

Bakta/Prokka emit no `gene` rows — only `CDS`, `tRNA`, `rRNA`, etc. `load_gff3` detects this automatically and treats each feature as gene-like, using `locus_tag` as the stable `gene_id`.

```python
db = GenomeDB()
db.load_genome("strain_BRG0069.fna")
db.load_gff3("strain_BRG0069.gff3", prefer="cds")   # cds is the right choice for prokaryotes

# locus_tag lookup
for g in db.get_genes_by_name("BRG0069_01234"):
    print(g.region(), g.product)
```

### 3. GenBank-only workflow

A single `.gbk` carries sequence + annotation. You don't need to call `load_genome` first:

```python
db = GenomeDB()
n = db.load_genbank("plasmid.gbk")       # auto-loads ORIGIN as genome, builds genes from CDS
print(f"{n} genes; {db.genome.get_genome_size():,} bp")

# Protein lookup
for pid, info in db.proteins.items():
    print(pid, "->", info["name"][:60], f"({len(info['sequence'])} aa)")
```

### 4. Overlap queries

```python
# "What overlaps this SNP?"
for g in db.genes_at("scaffold_1", 452_108):
    print(g.gene_id, g.region(), g.strand)

# "What's in this 10 kb window?"
for g in db.genes_overlapping("scaffold_1", 500_000, 510_000):
    print(g.as_tuple(with_strand=True))

# JSON-shaped summary for an API response
import json
print(json.dumps(db.describe_position("scaffold_1", 452_108), default=str, indent=2))
```

### 5. Export to a DataFrame

```python
import pandas as pd

arr, cols = db.genes_to_matrix(with_strand=True)
df = pd.DataFrame(arr, columns=cols)
print(df.head())
```

---

## GFF3 support notes

- **Hierarchical GFF3** (`gene` → `mRNA` → `exon`/`CDS`) is wired via `ID` and `Parent`. When a feature has multiple parents (comma-separated in `Parent=`), the block is attached to **each** parent transcript.
- **Gene-less GFF3** (Bakta/Prokka) — triggered when *no* row has `type == gene`. Each `CDS`/RNA row becomes a gene, using `locus_tag` (falling back to `ID` → `gene` → `Name`) as the stable `gene_id`, and the row itself is added as a single exon block under a pseudo-transcript so extraction still works.
- Attributes are URL-decoded (`%20` → space, etc.).
- InterProScan-style GFF3 (protein-coordinate features like `polypeptide`, `protein_match`) is detected and skipped with a warning — those coordinates are amino-acid, not genomic.
- **Tip:** bacterial genomes → `prefer="cds"`. Eukaryotic UTR-aware analyses → `prefer="exon"`.

---

## Performance

- Sequences stored as `np.uint8` ASCII bytes: **1 byte/base**. A 5 Mb bacterial genome is ~5 MB in RAM; a 3 Gb mammalian genome is ~3 GB.
- Region queries return zero-copy slices of the underlying array, then decode to `str` on output.
- Reverse-complement in `"bytes"` mode uses a 256-entry `uint8` lookup table (`_ASCII_RC`) + array reverse — no Python-level loops.
- Overlap queries: per-chromosome sorted `(starts, ends, gids)` arrays. `np.searchsorted` narrows to candidates that start ≤ `end`, then a vectorized mask filters by `ends >= start`. **O(log n)** candidate selection, **O(k)** filtering where `k` is the candidate count.

---

## Coordinate conventions

- 1-based, inclusive on every public method.
- `Genome.get_region(..., strand='-')` returns the reverse complement. Don't reverse-complement a second time at the call site.
- `get_region` **returns `None`** for an unknown chromosome and **raises `ValueError`** for an out-of-bounds region. Guard with `has_chromosome()` if that distinction matters.
- `genes_overlapping(chrom, start, end)` swaps arguments transparently if `start > end`.

---

## Project layout

```
GenomeManager/
  LICENSE                 # Apache-2.0
  NOTICE
  pyproject.toml
  requirements.txt
  README.md
  src/
    genomehandler/
      __init__.py         # exports: GenomeDB, Genome, Gene
      gene_class.py       # Gene (fields, exons/introns, extraction helpers)
      genome_class.py     # Genome (uint8 storage, get_region, reverse-complement)
      genome_db.py        # GenomeDB (GFF3/GenBank parsers, interval index, queries)
      genome_utils.py     # reverse_complement(str) helper
```

The `src/` layout keeps imports unambiguous during development and packaging.

---

## Troubleshooting

- **`Gene` object is not subscriptable.** Use attributes (`g.gene_id`) or `g.as_tuple()` — `Gene` is not a tuple/list.
- **Mixed class versions after edits.** Restart the kernel/process; cached `Gene` instances from a previous import will not have new attributes.
- **Out of bounds region.** Remember 1-based inclusive. Pass `strand` to `get_region` instead of reverse-complementing later.
- **`describe_position` / `get_genes_by_name` missing.** These require a version that includes the fix landed alongside this README. If you are pinned to an older revision, upgrade to `main` or ≥ v0.1.1.
- **Empty results from `get_genes_by_name`.** The multi-map is built from `Name`, `gene`, and `locus_tag` GFF3 attributes; if your file uses a non-standard attribute key, look the gene up by ID instead.
- **Loading a GBK after a FASTA.** `load_genbank(..., load_genome_if_present=True)` will **not** overwrite an already-loaded genome — it only fills the genome when none is loaded. Pass `False` to skip GBK sequence extraction entirely.

---

## Contributing

Issues and PRs are welcome at [ElytronBiotech/genomehandler](https://github.com/ElytronBiotech/genomehandler). For local dev:

```bash
git clone https://github.com/ElytronBiotech/genomehandler.git
cd genomehandler
python -m venv .venv && source .venv/bin/activate
pip install -e .
```

---

## License

© 2025 Elytron Biotech. Licensed under the **Apache License, Version 2.0** — see [LICENSE](./LICENSE) and [NOTICE](./NOTICE).
