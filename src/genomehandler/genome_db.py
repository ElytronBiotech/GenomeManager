from __future__ import annotations

import gzip
import os

import numpy as np
from Bio import SeqIO

from .genome_class import Genome
from .gene_class import Gene


class GenomeDB:
    """
    Minimal, fast, and clean genome database:
      - load_genome(): FASTA
      - load_gff3(): robust GFF3 (hierarchical gene/mRNA/exon/CDS or gene-less Bakta/Prokka)
      - genes_by_name: multi-map name -> [gene_id, ...]
      - get_genes_by_name(name): -> list[Gene]
      - genes_at(chrom, pos): -> list[Gene]
      - genes_overlapping(chrom, start, end): -> list[Gene]
      - load_genbank(): protein linkage (optional)
    """

    # ---------- construction ----------
    def __init__(self) -> None:
        self.genome: Genome = Genome()

        # core storage
        self.genes: dict[str, Gene] = {}                 # gene_id -> Gene
        self.genes_by_name: dict[str, list[str]] = {}    # name -> [gene_ids]

        # optional proteins
        self.proteins: dict[str, dict] = {}              # protein_id -> info
        self.protein_by_name: dict[str, dict] = {}       # product name -> info

        # fast interval index: chrom -> {"starts","ends","gids"} (np arrays)
        self._ivx: dict[str, dict[str, np.ndarray]] = {}

        # Track if genome has been loaded already
        self._genome_loaded: bool = False

    # ---------- helpers ----------
    @staticmethod
    def _open_any(path: str):
        """Open plain text or .gz transparently."""
        return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

    @staticmethod
    def _parse_gff3_attrs(attr_str: str) -> dict:
        """Parse `key=value;...` attributes with URL decoding."""
        from urllib.parse import unquote_plus
        out: dict[str, str] = {}
        if not attr_str or attr_str == ".":
            return out
        for field in attr_str.split(";"):
            if not field:
                continue
            if "=" in field:
                k, v = field.split("=", 1)
                out[k.strip()] = unquote_plus(v.strip())
            else:
                out[field.strip()] = "true"
        return out

    def _map_gene_name(self, name: str | None, gene_id: str) -> None:
        """Multi-map (one name → many gene_ids)."""
        if not name:
            return
        self.genes_by_name.setdefault(name, []).append(gene_id)

    def _build_interval_index(self) -> None:
        """Per-chromosome arrays for O(log n) overlap queries."""
        self._ivx.clear()
        by_chrom: dict[str, list[tuple[int, int, str]]] = {}
        for gid, g in self.genes.items():
            if not g.chromosome or g.start is None or g.end is None:
                continue
            by_chrom.setdefault(g.chromosome, []).append((g.start, g.end, gid))
        for chrom, ivals in by_chrom.items():
            ivals.sort(key=lambda x: (x[0], x[1]))
            starts = np.asarray([s for s, _, _ in ivals], dtype=np.int64)
            ends   = np.asarray([e for _, e, _ in ivals], dtype=np.int64)
            gids   = np.asarray([gid for *_, gid in ivals], dtype=object)
            self._ivx[chrom] = {"starts": starts, "ends": ends, "gids": gids}

    # ---------- public: genome ----------
    def load_genome(self, fasta_file: str) -> list[str]:
        """Load genome FASTA. Returns loaded chromosome names."""
        if not os.path.exists(fasta_file):
            raise FileNotFoundError(f"FASTA not found: {fasta_file}")
        names = self.genome.load_fasta(fasta_file)
        self._genome_loaded = True
        return names

    # ---------- public: GFF3 ----------
    def load_gff3(self, gff_path: str, *, prefer: str = "cds") -> int:
        """
        Load a GFF3 (supports 2 styles):
        A) "Rich": explicit 'gene' + 'mRNA' + 'exon'/'CDS'
        B) Bakta/Prokka "gene-less": only 'CDS'/RNA features w/ locus_tag

        prefer: "cds" | "exon"  (when both exon & CDS exist under transcripts)
        Returns number of genes loaded.
        """
        if prefer not in {"cds", "exon"}:
            raise ValueError("prefer must be 'cds' or 'exon'")
        if not os.path.exists(gff_path):
            raise FileNotFoundError(f"GFF3 not found: {gff_path}")

        # reset state that depends on annotation
        self.genes.clear()
        self.genes_by_name.clear()

        # quick sniff
        has_gene_features = False
        looks_like_interpro = False
        with self._open_any(gff_path) as sniff:
            for raw in sniff:
                if not raw:
                    continue
                if raw.startswith("##interproscan-version"):
                    looks_like_interpro = True
                if raw.startswith("#"):
                    continue
                cols = raw.rstrip("\n").split("\t", 8)
                if len(cols) != 9:
                    continue
                if cols[2].lower() == "gene":
                    has_gene_features = True
                    break

        # collectors for rich mode
        tx_to_gene: dict[str, str] = {}
        exons_by_tx: dict[str, list[tuple[int, int]]] = {}
        cds_by_tx: dict[str, list[tuple[int, int]]] = {}
        tx_strand: dict[str, str] = {}
        tx_chrom: dict[str, str] = {}
        pseudo_tx_for_gene: dict[str, str] = {}

        # which feature types we will treat as "gene-like" if there are no 'gene' rows
        gene_like = {
            "cds", "rrna", "trna", "ncrna", "srna", "tmrna", "misc_rna", "rna",
            "pseudogene"
        }
        # protein-centric rows to skip (coordinates are AA positions on proteins)
        protein_only = {"polypeptide", "protein", "protein_match"}

        # parse
        parsed_any_row = False
        with self._open_any(gff_path) as fh:
            for raw in fh:
                if not raw or raw.startswith("#"):
                    continue
                cols = raw.rstrip("\n").split("\t", 8)
                if len(cols) != 9:
                    continue

                seqid, source, ftype, start, end, score, strand, phase, attrs = cols
                t = ftype.lower()
                # Ignore InterProScan-style protein coordinates
                if t in protein_only:
                    continue
                print("$")
                try:
                    s_i = int(start)
                    e_i = int(end)
                except ValueError:
                    continue

                a = self._parse_gff3_attrs(attrs)
                chrom = seqid  # keep downstream naming consistent

                parsed_any_row = True

                if has_gene_features:
                    # ---- A) gene/mRNA/exon/CDS wiring ----
                    if t == "gene":
                        gid = a.get("ID") or a.get("gene_id") or f"gene:{chrom}:{s_i}-{e_i}:{strand}"
                        gname = a.get("Name") or a.get("gene_name") or a.get("locus_tag") or a.get("gene") or gid
                        if gid not in self.genes:
                            print("#")
                            print(gname)    
                            g = Gene(gid, gname, chrom, s_i, e_i, strand)
                            self.genes[gid] = g
                        else:
                            g = self.genes[gid]
                            if not g.chromosome:
                                g.chromosome = chrom
                            if g.start is None:
                                g.start = s_i
                            if g.end is None:
                                g.end = e_i
                            if not g.strand:
                                g.strand = strand
                        
                        self._map_gene_name(gname, gid)
                        if a.get("locus_tag"):
                            self._map_gene_name(a["locus_tag"], gid)
                        # keep collecting; no 'continue' needed

                    elif t in {"mrna", "transcript"}:
                        tid = a.get("ID")
                        parent = (a.get("Parent") or "").split(",")[0]
                        if tid and parent:
                            tx_to_gene[tid] = parent
                            tx_strand[tid] = strand
                            tx_chrom[tid] = chrom
                            exons_by_tx.setdefault(tid, [])
                            cds_by_tx.setdefault(tid, [])

                    elif t in {"exon", "cds"}:
                        parents = (a.get("Parent") or "").split(",")
                        for pid in parents:
                            if not pid:
                                continue
                            if pid in tx_to_gene:
                                if t == "exon":
                                    exons_by_tx.setdefault(pid, []).append((s_i, e_i))
                                else:
                                    cds_by_tx.setdefault(pid, []).append((s_i, e_i))
                                tx_strand[pid] = strand
                                tx_chrom[pid] = chrom
                            elif pid in self.genes:
                                # feature attached directly to gene (no transcript feature)
                                tid = pseudo_tx_for_gene.setdefault(pid, f"{pid}.tx")
                                tx_to_gene.setdefault(tid, pid)
                                tx_strand.setdefault(tid, strand)
                                tx_chrom.setdefault(tid, chrom)
                                if t == "exon":
                                    exons_by_tx.setdefault(tid, []).append((s_i, e_i))
                                else:
                                    cds_by_tx.setdefault(tid, []).append((s_i, e_i))

                else:
                    # ---- B) gene-less (Bakta/Prokka) ----
                    if t not in gene_like:
                        continue

                    gid = (
                        a.get("locus_tag")
                        or a.get("ID")
                        or a.get("gene")
                        or a.get("Name")
                        or f"feat:{chrom}:{s_i}-{e_i}:{strand}"
                    )
                    gname = a.get("Name") or a.get("gene_name") or a.get("locus_tag") or a.get("gene") or gid
                    print("%")
                    print(gname)
                    g = self.genes.get(gid)
                    if not g:
                        g = Gene(gid, gname, chrom, s_i, e_i, strand)
                        self.genes[gid] = g
                        self._map_gene_name(gname, gid)
                        if a.get("locus_tag"):
                            self._map_gene_name(a["locus_tag"], gid)
                    else:
                        # grow bounds if same locus spans multiple parts
                        if g.start is None or s_i < g.start:
                            g.start = s_i
                        if g.end is None or e_i > g.end:
                            g.end = e_i
                        if not g.chromosome:
                            g.chromosome = chrom
                        if not g.strand:
                            g.strand = strand

                    # treat each piece as an exon-like block (enables downstream extraction)
                    tid = pseudo_tx_for_gene.setdefault(gid, f"{gid}.tx")
                    g.add_exon(s_i, e_i, tid)
                    if hasattr(g, "product") and a.get("product"):
                        g.product = a["product"]

        # finalize transcripts (rich mode)
        if has_gene_features and self.genes:
            for tid, gid in tx_to_gene.items():
                gene = self.genes.get(gid)
                if not gene:
                    continue
                blocks = exons_by_tx.get(tid, [])
                if prefer == "cds" and cds_by_tx.get(tid):
                    blocks = cds_by_tx[tid]
                if not blocks:
                    continue
                for s, e in sorted(blocks, key=lambda x: (x[0], x[1])):
                    gene.add_exon(s, e, tid)

        # per-gene post-processing if Gene implements helpers
        for g in self.genes.values():
            if hasattr(g, "set_gene_start_end"):
                g.set_gene_start_end()
            if hasattr(g, "calculate_introns"):
                g.calculate_introns()

        # build interval index keyed by chrom
        self._build_interval_index()

        # friendly warnings
        if not self.genes and looks_like_interpro:
            print("Warning: GFF3 looks like InterProScan (protein-centric). Skipping protein features; no genomic genes loaded.")
        elif not parsed_any_row:
            print("Warning: No parseable rows found in GFF3.")

        return len(self.genes)


    def get_genes_by_id(self, gene_id: str) -> list["Gene"]:
        """Return all Gene objects that match this gene_id (exact key)."""
        g = self.genes.get(gene_id)
        return [g] if g else []


    def genes_overlapping(self, chrom: str, start: int, end: int) -> list[Gene]:
        """Return ALL genes overlapping the closed interval [start, end]."""
        if start > end:
            start, end = end, start
        iv = self._ivx.get(chrom)
        if not iv:
            return []
        starts, ends, gids = iv["starts"], iv["ends"], iv["gids"]
        k = np.searchsorted(starts, end, side="right")
        if k == 0:
            return []
        mask = ends[:k] >= start
        if not np.any(mask):
            return []
        return [self.genes[gids[i]] for i in np.nonzero(mask)[0]]

    def describe_position(self, chrom: str, pos: int) -> list[dict]:
        """
        Convenience: describe all genes overlapping a position.
        Returns a list of dicts (one per gene). Empty list if none.
        """
        items: list[dict] = []
        for g in self.genes_at(chrom, pos):
            # compute cDNA on-demand if your Gene supports it and genome is loaded
            if (getattr(g, "cdna_sequence", None) is None
                and self.genome.get_chromosome(g.chromosome) is not None):
                if hasattr(g, "extract_cdna_sequence"):
                    g.extract_cdna_sequence(self.genome)
            items.append({
                "gene_id": g.gene_id,
                "gene_name": g.gene_name,
                "location": f"{g.chromosome}:{g.start}-{g.end}",
                "start": g.start,
                "end": g.end,
                "strand": g.strand,
                "exons": getattr(g, "exons", None),
                "introns": getattr(g, "introns", None),
                "cdna": getattr(g, "cdna_sequence", None),
                "protein_id": getattr(g, "protein_id", None),
            })
        return items

    def _load_genome_from_seqs(self, seqs):
        if not seqs:
            return
        self.genome.load_sequences(seqs)  # no temp file
        self._genome_loaded = True

    # ---------- optional: genbank proteins ----------
    def load_genbank(self, gbk_file: str, *, load_genome_if_present: bool = True) -> int:
        """
        Load CDS features from GenBank to collect protein sequences and link them to genes.
        If self.genes is empty, also build genes from GBK ('gene' and/or 'CDS' features).

        If load_genome_if_present is True and the GBK contains sequence records (ORIGIN),
        the genome will be auto-loaded (only if not already loaded).
        """
        import os, gzip
        from Bio import SeqIO
        from Bio.SeqFeature import CompoundLocation
        print("Loading GenBank:", gbk_file)
        if not os.path.exists(gbk_file):
            raise FileNotFoundError(f"GenBank not found: {gbk_file}")

        # Clear protein indices; (re)build from GBK below
        self.proteins.clear()
        self.protein_by_name.clear()

        open_fn = gzip.open if gbk_file.endswith(".gz") else open
        with open_fn(gbk_file, "rt") as handle:
            records = list(SeqIO.parse(handle, "genbank"))

        # Optional genome autoload from GBK sequence, if not already loaded
        if load_genome_if_present and not getattr(self, "_genome_loaded", False):
            seqs = []
            for rec in records:
                if getattr(rec, "seq", None) and len(rec.seq) > 0:
                    seqs.append((rec.id or rec.name or f"record_{len(seqs)+1}",
                                str(rec.seq).upper()))
            if seqs:
                self._load_genome_from_seqs(seqs)  # should set self._genome_loaded internally

        # Helper: convert Biopython location to 1-based inclusive coords + strand char
        def _loc_info(loc):
            s = int(loc.start) + 1   # Biopython is 0-based, GFF/GBK locations are 1-based
            e = int(loc.end)         # inclusive end
            strand = "+" if getattr(loc, "strand", 0) == 1 else ("-" if getattr(loc, "strand", 0) == -1 else ".")
            return s, e, strand

        # If genes are missing, build them from GBK
        need_build_genes = (len(self.genes) == 0)

        if need_build_genes:
            # Pass 1: create genes from explicit 'gene' features (best source of names & bounds)
            for rec in records:
                chrom = rec.id or rec.name or "chr"
                for feat in getattr(rec, "features", []):
                    if feat.type != "gene":
                        continue
                    q = feat.qualifiers
                    gid = q.get("locus_tag", [""])[0]
                    gname = q.get("gene", [""])[0] or gid
                    s_i, e_i, strand = _loc_info(feat.location)
                    if gid not in self.genes:
                        g = Gene(gid, gname, chrom, s_i, e_i, strand)
                        self.genes[gid] = g
                    else:
                        g = self.genes[gid]
                        # be conservative updating if something is missing
                        if not getattr(g, "chromosome", None):
                            g.chromosome = chrom
                        if getattr(g, "start", None) is None: g.start = s_i
                        if getattr(g, "end", None)   is None: g.end   = e_i
                        if not getattr(g, "strand", None):    g.strand = strand

                    # name indices
                    if gname:
                        self._map_gene_name(gname, gid)

            # Pass 2: ensure any CDS without a preceding 'gene' still creates a gene,
            # and add exon-like blocks for each CDS (handles joined locations)
            for rec in records:
                chrom = rec.id or rec.name or "chr"
                for feat in getattr(rec, "features", []):
                    if feat.type != "CDS":
                        continue
                    q = feat.qualifiers
                    locus = q.get("locus_tag", [""])[0]
                    gsym  = q.get("gene", [""])[0]
                    gid   = locus or gsym or q.get("old_locus_tag", [""])[0] or f"feat:{chrom}:{feat.location.start}-{feat.location.end}"
                    gname = gsym or locus or gid

                    # create missing gene from the CDS bounds
                    if gid not in self.genes:
                        s_i, e_i, strand = _loc_info(feat.location)
                        g = Gene(gid, gname, chrom, s_i, e_i, strand)
                        self.genes[gid] = g
                        if gname:
                            self._map_gene_name(gname, gid)
                        if locus:
                            self._map_gene_name(locus, gid)
                    else:
                        g = self.genes[gid]

                    # add exon-like blocks for this CDS
                    if isinstance(feat.location, CompoundLocation):
                        for part in feat.location.parts:
                            ps, pe, _ = _loc_info(part)
                            g.add_exon(ps, pe, f"{gid}.tx")
                    else:
                        s_i, e_i, _ = _loc_info(feat.location)
                        g.add_exon(s_i, e_i, f"{gid}.tx")

                    # stash a useful product if present
                    prod = q.get("product", [""])[0]
                    if prod and hasattr(g, "product") and not getattr(g, "product", None):
                        g.product = prod

        # Always attach protein info from CDS features
        for rec in records:
            for feat in getattr(rec, "features", []):
                if feat.type != "CDS":
                    continue
                q = feat.qualifiers
                protein_id   = q.get("protein_id", [""])[0]
                protein_name = q.get("product", [""])[0]
                gene_symbol  = q.get("gene", [""])[0]
                protein_seq  = q.get("translation", [""])[0]

                if protein_id:
                    info = {
                        "id": protein_id,
                        "name": protein_name,
                        "gene_name": gene_symbol,
                        "sequence": protein_seq,
                    }
                    self.proteins[protein_id] = info
                    if protein_name:
                        self.protein_by_name[protein_name] = info

                # attach to a matching Gene if we can
                candidates = []
                locus = q.get("locus_tag", [""])[0]
                if locus:
                    # Prefer an exact locus_tag match
                    g = self.genes.get(locus)
                    if g:
                        candidates = [g]
                if not candidates and gene_symbol:
                    for gid in self.genes_by_name.get(gene_symbol, []):
                        if gid in self.genes:
                            candidates.append(self.genes[gid])

                # If still nothing and we just built genes, fall back to making one quickly
                if not candidates and need_build_genes:
                    chrom = rec.id or rec.name or "chr"
                    s_i, e_i, strand = _loc_info(feat.location)
                    gid = locus or gene_symbol or f"feat:{chrom}:{feat.location.start}-{feat.location.end}"
                    gname = gene_symbol or locus or gid
                    g = self.genes.get(gid)
                    if not g:
                        g = Gene(gid, gname, chrom, s_i, e_i, strand)
                        self.genes[gid] = g
                        if gname:
                            self._map_gene_name(gname, gid)
                        if locus:
                            self._map_gene_name(locus, gid)
                    candidates = [g]

                for g in candidates:
                    if protein_seq:
                        g.protein_sequence = protein_seq
                    if protein_id:
                        g.protein_id = protein_id
                    # ensure at least one CDS block exists
                    if isinstance(feat.location, CompoundLocation):
                        for part in feat.location.parts:
                            ps, pe, _ = _loc_info(part)
                            g.add_exon(ps, pe, f"{g.gene_id}.tx" if hasattr(g, "gene_id") else f"{g.id}.tx")
                    else:
                        s_i, e_i, _ = _loc_info(feat.location)
                        g.add_exon(s_i, e_i, f"{g.gene_id}.tx" if hasattr(g, "gene_id") else f"{g.id}.tx")

        # Post-process genes if helpers exist
        for g in self.genes.values():
            if hasattr(g, "set_gene_start_end"):
                g.set_gene_start_end()
            if hasattr(g, "calculate_introns"):
                g.calculate_introns()

        # Refresh interval index after potential additions
        self._build_interval_index()

        # Return number of genes now known
        return len(self.genes)
