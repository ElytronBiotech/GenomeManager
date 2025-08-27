# GenomeHandler

Minimal, fast, and clean Python library for working with genomes and gene annotations.

- **Inputs:** FASTA (optionally `.gz`), GFF3 (hierarchical or Bakta/Prokka-style “gene-less”), optional GenBank for protein info.
- **Core classes:** `GenomeDB` (queries + I/O), `Genome` (memory‑efficient sequence store), `Gene` (structure + convenience helpers).
- **Coordinates:** 1‑based, inclusive (e.g., `start=1, end=3` → 3 bases).
- **Performance:** 1 byte/base storage with NumPy, vectorized interval index for O(log n) overlap queries.

> Requires Python **3.11+**, NumPy, Biopython.

---

## Installation

### From GitHub (latest `main`)
```bash
pip install "genomehandler @ git+https://github.com/ElytronBiotech/genomehandler.git@main"
```

### From a tagged release
```bash
pip install "genomehandler @ git+https://github.com/ElytronBiotech/genomehandler.git@v0.1.0"
```

### Local (editable) install during development
```bash
python -m venv .venv && source .venv/bin/activate    # or use conda/mamba
pip install -U pip
pip install -e .
```

---

## Quickstart

```python
from genomehandler import GenomeDB

db = GenomeDB()

# 1) Load genome (FASTA or FASTA.gz)
db.load_genome("R000044.fasta.gz")

# 2) Load annotations (GFF3; works with/without explicit 'gene')
db.load_gff3("R000044.gff3", prefer="cds")  # or prefer="exon"

# 3) Look up by name (multi‑mapping): returns List[Gene]
genes = db.get_genes_by_name("tmk")
for g in genes:
    print(g.as_tuple(with_strand=True), g.region())

# 4) Overlap queries
hits = db.genes_at("scaffold_1", 5600)                 # List[Gene]
span = db.genes_overlapping("scaffold_1", 5200, 10500) # List[Gene]

# 5) cDNA extraction (strand‑aware slicing under the hood)
g = genes[0]
seq = g.extract_cdna_sequence(db.genome)
print(g.region(), len(seq))
```

---

## API Overview

### `GenomeDB`
- `load_genome(fasta: str) -> list[str]` — load chromosomes from FASTA/FASTA.gz.
- `load_gff3(path: str, *, prefer: str = "cds") -> int` — load GFF3 (hierarchy or gene‑less). `prefer` switches exon vs CDS when both exist.
- `get_genes_by_name(name: str) -> list[Gene]` — name/symbol/locus_tag → Gene objects.
- `genes_at(chrom: str, pos: int) -> list[Gene]` — all genes covering a position.
- `genes_overlapping(chrom: str, start: int, end: int) -> list[Gene]` — all genes intersecting a span.
- `describe_position(chrom: str, pos: int) -> list[dict]` — JSON‑ready convenience summary.
- `load_genbank(path: str) -> int` — (optional) add protein info and link by gene symbol.

### `Genome`
- `load_fasta(path: str) -> list[str]` — supports `.gz`.
- `get_region(chrom: str, start: int, end: int, *, strand: str = "+") -> str | None` — 1‑based inclusive; reverse‑complement when `strand='-'`.
- `get_chromosome_names()`, `get_chromosome_size(name)`, `get_genome_size()`.

### `Gene`
- Fields: `gene_id`, `gene_name`, `chromosome`, `start`, `end`, `strand`, `transcripts`, `exons`, `introns`.
- Methods: `add_exon(...)`, `calculate_introns()`, `set_gene_start_end()`, `extract_cdna_sequence(genome, transcript=None)`.
- **Convenience:** `as_tuple(with_strand=False)`, `region()`, `exons_as_tuples()`, `introns_as_tuples()`.

---

## Coordinate Conventions

- All public methods use **1‑based, inclusive** coordinates.
- Negative‑strand sequences are returned reverse‑complemented by `Genome.get_region(..., strand='-')`.
- GFF3 parsing preserves strand and concatenates blocks in genomic order.

---

## GFF3 Support Notes

- **Hierarchical** files with `gene`/`mRNA`/`exon`/`CDS` are wired via `Parent`/`ID`.
- **Gene‑less** (Bakta/Prokka) files are supported by treating `CDS`/RNA rows as gene‑like features:
  - stable `gene_id` from `locus_tag` (or `ID`, `gene`, `Name`),
  - a single pseudo‑transcript per gene,
  - each segment added as a block (`add_exon`) so extraction works uniformly.
- Attributes are URL‑decoded; multi‑parent exons are attached to each transcript.

**Tip:** For bacterial genomes use `prefer="cds"` (coding sequence). For cDNA/UTR‑aware analyses use `prefer="exon"`.

---

## Performance

- FASTA stored as ASCII bytes via NumPy: **1 byte/base**, fast views for slicing.
- Overlap queries use per‑chromosome interval arrays (`starts`, `ends`, `ids`) plus `np.searchsorted` for O(log n) candidate selection.

---

## Project Layout

```
genomehandler/
  docs/
  src/
    genomehandler/
      __init__.py
      gene_class.py
      genome_class.py
      genome_db.py
      genome_utils.py
```

> The `src/` layout avoids import path confusion during development and packaging.

---

## Requirements

- Python **3.11+**
- `numpy>=1.24`
- `biopython>=1.81`

Install dev tools with extras or requirements files if provided.

---

## Troubleshooting

- **“Gene object is not subscriptable.”** Use attributes or `g.as_tuple()` instead of `g[0]`.
- **Mixed class versions / missing attributes.** Restart the kernel/process and reload modules so all `Gene` objects use the current class.
- **Out of bounds region.** Remember 1‑based inclusive coordinates; pass `strand` to `get_region` instead of reverse‑complementing later.

---

## License

© 2025 Elytron Biotech. Licensed under **Apache‑2.0**.