"""Tests for the Gene class (structure, cDNA extraction)."""

import pytest

from genomehandler import Gene, Genome


@pytest.fixture
def genome():
    """Genome with a known sequence for predictable cDNA tests."""
    g = Genome()
    #        positions: 1234567890123456789012345678901234567890
    g.load_sequences([("chr1", "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCC")])
    return g


class TestGeneConstruction:
    def test_positional_args(self):
        g = Gene("g1", "geneA", "chr1", 100, 200, "+")
        assert g.gene_id == "g1"
        assert g.gene_name == "geneA"
        assert g.chromosome == "chr1"
        assert g.start == 100
        assert g.end == 200
        assert g.strand == "+"

    def test_defaults(self):
        g = Gene("g1", "geneA", "chr1")
        assert g.start is None
        assert g.end is None
        assert g.strand is None
        assert g.exons == {}
        assert g.transcripts == []


class TestAddExon:
    def test_add_exon_default_isoform(self):
        g = Gene("g1", "geneA", "chr1", 1, 100, "+")
        g.add_exon(10, 20)
        assert (10, 20) in g.exons
        assert g.exons[(10, 20)] == ["g1"]
        assert "g1" in g.transcripts

    def test_add_exon_named_isoform(self):
        g = Gene("g1", "geneA", "chr1", 1, 100, "+")
        g.add_exon(10, 20, "tx1")
        g.add_exon(30, 40, "tx1")
        assert g.transcripts == ["tx1"]
        assert len(g.exons) == 2

    def test_duplicate_exon_isoform(self):
        g = Gene("g1", "geneA", "chr1", 1, 100, "+")
        g.add_exon(10, 20, "tx1")
        g.add_exon(10, 20, "tx1")
        assert g.exons[(10, 20)] == ["tx1"]


class TestCalculateIntrons:
    def test_single_intron(self):
        g = Gene("g1", "geneA", "chr1", 1, 100, "+")
        g.add_exon(10, 20, "tx1")
        g.add_exon(31, 40, "tx1")
        introns = g.calculate_introns()
        assert (21, 30) in introns
        assert "tx1" in introns[(21, 30)]

    def test_no_introns_adjacent(self):
        g = Gene("g1", "geneA", "chr1", 1, 100, "+")
        g.add_exon(10, 20, "tx1")
        g.add_exon(21, 30, "tx1")
        introns = g.calculate_introns()
        assert len(introns) == 0


class TestSetGeneStartEnd:
    def test_set_from_exons(self):
        g = Gene("g1", "geneA", "chr1", None, None, "+")
        g.add_exon(50, 60, "tx1")
        g.add_exon(80, 90, "tx1")
        g.set_gene_start_end()
        assert g.start == 50
        assert g.end == 90

    def test_no_exons(self):
        g = Gene("g1", "geneA", "chr1", 1, 100, "+")
        result = g.set_gene_start_end()
        assert result is None


class TestExtractCdna:
    def test_plus_strand_single_exon(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 4, "+")
        g.add_exon(1, 4)
        cdna = g.extract_cdna_sequence(genome)
        assert cdna == "AAAA"

    def test_plus_strand_multi_exon(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 12, "+")
        g.add_exon(1, 4)   # AAAA
        g.add_exon(9, 12)  # GGGG
        cdna = g.extract_cdna_sequence(genome)
        assert cdna == "AAAAGGGG"

    def test_minus_strand_single_exon(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 4, "-")
        g.add_exon(1, 4)
        cdna = g.extract_cdna_sequence(genome)
        # RC of AAAA = TTTT
        assert cdna == "TTTT"

    def test_minus_strand_multi_exon_order(self, genome):
        """Regression test for bug 1a: multi-exon minus-strand order."""
        g = Gene("g1", "geneA", "chr1", 1, 12, "-")
        g.add_exon(1, 4)   # AAAA
        g.add_exon(9, 12)  # GGGG
        cdna = g.extract_cdna_sequence(genome)
        # Correct: join AAAA + GGGG = AAAAGGGG, then RC whole = CCCCTTTT
        assert cdna == "CCCCTTTT"

    def test_transcript_specific(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 12, "+")
        g.add_exon(1, 4, "tx1")
        g.add_exon(5, 8, "tx2")
        g.add_exon(9, 12, "tx1")
        cdna = g.extract_cdna_sequence(genome, transcript="tx1")
        assert cdna == "AAAAGGGG"  # exons 1-4 + 9-12

    def test_unknown_transcript_returns_none(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 4, "+")
        g.add_exon(1, 4)
        assert g.extract_cdna_sequence(genome, transcript="nope") is None

    def test_no_exons_returns_none(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 4, "+")
        assert g.extract_cdna_sequence(genome) is None

    def test_caches_when_no_transcript(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 4, "+")
        g.add_exon(1, 4)
        g.extract_cdna_sequence(genome)
        assert g.cdna_sequence == "AAAA"

    def test_does_not_cache_transcript_specific(self, genome):
        g = Gene("g1", "geneA", "chr1", 1, 8, "+")
        g.add_exon(1, 4, "tx1")
        g.add_exon(5, 8, "tx1")
        g.extract_cdna_sequence(genome, transcript="tx1")
        assert g.cdna_sequence is None


class TestHelpers:
    def test_as_tuple(self):
        g = Gene("g1", "geneA", "chr1", 10, 20, "+")
        assert g.as_tuple() == ("g1", "geneA", "chr1", 10, 20)
        assert g.as_tuple(with_strand=True) == ("g1", "geneA", "chr1", 10, 20, "+")

    def test_region(self):
        g = Gene("g1", "geneA", "chr1", 10, 20, "+")
        assert g.region() == "chr1:10-20"

    def test_region_unknown(self):
        g = Gene("g1", "geneA", "chr1")
        assert g.region() == "chr1:<unknown>"

    def test_str(self):
        g = Gene("g1", "geneA", "chr1", 10, 20, "+")
        s = str(g)
        assert "geneA" in s
        assert "g1" in s

    def test_call(self):
        g = Gene("g1", "geneA", "chr1", 10, 20, "+")
        assert g() == g.as_tuple()
