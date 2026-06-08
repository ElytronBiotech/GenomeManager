"""Tests for GenBank loading in GenomeDB."""

import pytest

from genomehandler import GenomeDB


def _write_genbank(path, records_text):
    """Write a GenBank file from raw text."""
    with open(path, "w") as f:
        f.write(records_text)


SIMPLE_GBK = """\
LOCUS       contig1                   40 bp    DNA     linear   BCT 01-JAN-2025
DEFINITION  Test contig.
ACCESSION   contig1
VERSION     contig1.1
FEATURES             Location/Qualifiers
     gene            5..20
                     /locus_tag="LOC_001"
                     /gene="testGene"
     CDS             5..20
                     /locus_tag="LOC_001"
                     /gene="testGene"
                     /protein_id="WP_000001.1"
                     /product="test protein"
                     /translation="MFGH"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""

COMPOUND_GBK = """\
LOCUS       contig1                   40 bp    DNA     linear   BCT 01-JAN-2025
DEFINITION  Test contig with compound location.
ACCESSION   contig1
VERSION     contig1.1
FEATURES             Location/Qualifiers
     gene            complement(join(5..12,21..28))
                     /locus_tag="LOC_002"
                     /gene="joinGene"
     CDS             complement(join(5..12,21..28))
                     /locus_tag="LOC_002"
                     /gene="joinGene"
                     /protein_id="WP_000002.1"
                     /product="joined protein"
                     /translation="MKRL"
ORIGIN
        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt
//
"""


class TestLoadGenbank:
    def test_basic_load(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        n = db.load_genbank(str(gbk))
        assert n >= 1
        assert "LOC_001" in db.genes

    def test_genome_autoload(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        # Biopython uses VERSION (rec.id) as the chromosome name
        assert db.genome.has_chromosome("contig1.1")
        assert db.genome.get_chromosome_size("contig1.1") == 40

    def test_genome_autoload_disabled(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk), load_genome_if_present=False)
        assert not db.genome.has_chromosome("contig1")

    def test_protein_linkage(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        g = db.genes["LOC_001"]
        assert g.protein_id == "WP_000001.1"
        assert g.protein_sequence == "MFGH"
        assert g.product == "test protein"

    def test_protein_table(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        assert "WP_000001.1" in db.proteins
        info = db.proteins["WP_000001.1"]
        assert info["name"] == "test protein"
        assert info["sequence"] == "MFGH"

    def test_compound_location_exons(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), COMPOUND_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        g = db.genes["LOC_002"]
        exon_coords = g.exons_as_tuples()
        assert (5, 12) in exon_coords
        assert (21, 28) in exon_coords

    def test_compound_minus_strand(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), COMPOUND_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        g = db.genes["LOC_002"]
        assert g.strand == "-"

    def test_name_lookup(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        genes = db.get_genes_by_name("testGene")
        assert len(genes) >= 1
        assert genes[0].gene_id == "LOC_001"

    def test_file_not_found(self):
        db = GenomeDB()
        with pytest.raises(FileNotFoundError):
            db.load_genbank("/nonexistent/path.gbk")

    def test_cdna_extraction_after_genbank(self, tmp_path):
        gbk = tmp_path / "test.gbk"
        _write_genbank(str(gbk), SIMPLE_GBK)

        db = GenomeDB()
        db.load_genbank(str(gbk))
        g = db.genes["LOC_001"]
        cdna = g.extract_cdna_sequence(db.genome)
        assert cdna is not None
        assert len(cdna) == 16  # positions 5-20 inclusive
