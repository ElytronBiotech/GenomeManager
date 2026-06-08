"""Tests for GFF3 loading in GenomeDB."""

import pytest

from genomehandler import GenomeDB


def _write_fasta(path, seqs):
    """Write a simple FASTA file."""
    with open(path, "w") as f:
        for name, seq in seqs:
            f.write(f">{name}\n{seq}\n")


def _write_gff3(path, lines):
    """Write a GFF3 file from a list of tab-separated lines."""
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for line in lines:
            f.write(line + "\n")


class TestHierarchicalGff3:
    @pytest.fixture
    def db_with_hierarchical(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("chr1", "A" * 1000)])

        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), [
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=geneA",
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1",
            "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1",
            "chr1\t.\texon\t300\t400\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t110\t190\t.\t+\t0\tParent=tx1",
            "chr1\t.\tCDS\t310\t390\t.\t+\t0\tParent=tx1",
        ])

        db = GenomeDB()
        db.load_genome(str(fasta))
        db.load_gff3(str(gff), prefer="exon")
        return db

    def test_loads_gene(self, db_with_hierarchical):
        db = db_with_hierarchical
        assert "gene1" in db.genes
        g = db.genes["gene1"]
        assert g.gene_name == "geneA"
        assert g.strand == "+"

    def test_prefer_exon(self, db_with_hierarchical):
        db = db_with_hierarchical
        g = db.genes["gene1"]
        exon_coords = g.exons_as_tuples()
        assert (100, 200) in exon_coords
        assert (300, 400) in exon_coords

    def test_prefer_cds(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("chr1", "A" * 1000)])

        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), [
            "chr1\t.\tgene\t100\t500\t.\t+\t.\tID=gene1;Name=geneA",
            "chr1\t.\tmRNA\t100\t500\t.\t+\t.\tID=tx1;Parent=gene1",
            "chr1\t.\texon\t100\t200\t.\t+\t.\tParent=tx1",
            "chr1\t.\texon\t300\t400\t.\t+\t.\tParent=tx1",
            "chr1\t.\tCDS\t110\t190\t.\t+\t0\tParent=tx1",
            "chr1\t.\tCDS\t310\t390\t.\t+\t0\tParent=tx1",
        ])

        db = GenomeDB()
        db.load_genome(str(fasta))
        db.load_gff3(str(gff), prefer="cds")
        g = db.genes["gene1"]
        exon_coords = g.exons_as_tuples()
        assert (110, 190) in exon_coords
        assert (310, 390) in exon_coords

    def test_name_lookup(self, db_with_hierarchical):
        db = db_with_hierarchical
        genes = db.get_genes_by_name("geneA")
        assert len(genes) == 1
        assert genes[0].gene_id == "gene1"

    def test_unknown_name(self, db_with_hierarchical):
        db = db_with_hierarchical
        assert db.get_genes_by_name("nonexistent") == []


class TestPreferExonFallback:
    """Regression tests for bug 1b: prefer='exon' with only CDS blocks."""

    def test_fallback_to_cds_when_no_exon_features(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("chr1", "ACGT" * 250)])

        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), [
            "chr1\t.\tgene\t100\t400\t.\t+\t.\tID=gene1;Name=geneA",
            "chr1\t.\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene1",
            "chr1\t.\tCDS\t110\t190\t.\t+\t0\tParent=tx1",
            "chr1\t.\tCDS\t310\t390\t.\t+\t0\tParent=tx1",
        ])

        db = GenomeDB()
        db.load_genome(str(fasta))
        db.load_gff3(str(gff), prefer="exon")
        g = db.genes["gene1"]
        # Should fall back to CDS blocks since no exon features exist
        assert len(g.exons) > 0
        exon_coords = g.exons_as_tuples()
        assert (110, 190) in exon_coords
        assert (310, 390) in exon_coords


class TestGeneLessGff3:
    """Tests for Bakta/Prokka-style gene-less GFF3."""

    @pytest.fixture
    def db_geneless(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("contig1", "ACGTACGT" * 100)])

        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), [
            "contig1\t.\tCDS\t10\t100\t.\t+\t0\tID=cds1;locus_tag=LOCUS_001;product=hypothetical",
            "contig1\t.\tCDS\t200\t300\t.\t-\t0\tID=cds2;locus_tag=LOCUS_002;gene=geneB",
        ])

        db = GenomeDB()
        db.load_genome(str(fasta))
        db.load_gff3(str(gff))
        return db

    def test_genes_created(self, db_geneless):
        db = db_geneless
        assert "LOCUS_001" in db.genes
        assert "LOCUS_002" in db.genes

    def test_locus_tag_lookup(self, db_geneless):
        db = db_geneless
        genes = db.get_genes_by_name("LOCUS_001")
        assert len(genes) >= 1

    def test_gene_symbol_lookup(self, db_geneless):
        db = db_geneless
        # In gene-less mode, gname = Name || gene_name || locus_tag || gene
        # With only locus_tag and gene attrs, locus_tag wins for name mapping
        genes = db.get_genes_by_name("LOCUS_002")
        assert len(genes) >= 1
        assert genes[0].gene_id == "LOCUS_002"

    def test_product_set(self, db_geneless):
        db = db_geneless
        g = db.genes["LOCUS_001"]
        assert g.product == "hypothetical"


class TestIntervalQueries:
    @pytest.fixture
    def db(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("chr1", "A" * 1000)])

        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), [
            "chr1\t.\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=g1",
            "chr1\t.\tgene\t150\t300\t.\t-\t.\tID=gene2;Name=g2",
            "chr1\t.\tgene\t500\t600\t.\t+\t.\tID=gene3;Name=g3",
        ])

        db = GenomeDB()
        db.load_genome(str(fasta))
        db.load_gff3(str(gff))
        return db

    def test_genes_at_overlap(self, db):
        hits = db.genes_at("chr1", 175)
        ids = {g.gene_id for g in hits}
        assert "gene1" in ids
        assert "gene2" in ids

    def test_genes_at_no_hit(self, db):
        assert db.genes_at("chr1", 400) == []

    def test_genes_overlapping_range(self, db):
        hits = db.genes_overlapping("chr1", 250, 550)
        ids = {g.gene_id for g in hits}
        assert "gene2" in ids
        assert "gene3" in ids

    def test_genes_overlapping_swaps_args(self, db):
        # start > end should still work
        hits = db.genes_overlapping("chr1", 550, 250)
        ids = {g.gene_id for g in hits}
        assert "gene2" in ids

    def test_unknown_chrom_empty(self, db):
        assert db.genes_at("chrX", 100) == []

    def test_invalid_prefer_raises(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("chr1", "AAAA")])
        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), ["chr1\t.\tgene\t1\t4\t.\t+\t.\tID=g1"])

        db = GenomeDB()
        db.load_genome(str(fasta))
        with pytest.raises(ValueError):
            db.load_gff3(str(gff), prefer="invalid")


class TestDescribePosition:
    def test_describe_with_genome(self, tmp_path):
        fasta = tmp_path / "genome.fasta"
        _write_fasta(str(fasta), [("chr1", "ACGT" * 50)])

        gff = tmp_path / "ann.gff3"
        _write_gff3(str(gff), [
            "chr1\t.\tCDS\t5\t20\t.\t+\t0\tID=cds1;locus_tag=LOC1",
        ])

        db = GenomeDB()
        db.load_genome(str(fasta))
        db.load_gff3(str(gff))

        items = db.describe_position("chr1", 10)
        assert len(items) == 1
        assert items[0]["gene_id"] == "LOC1"
        assert items[0]["cdna"] is not None
