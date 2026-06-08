"""Tests for the Genome class (sequence storage and retrieval)."""

import pytest

from genomehandler import Genome


@pytest.fixture
def genome():
    g = Genome()
    # chr1: 20 bases, mixed case to verify case preservation
    g.load_sequences([("chr1", "ACGTACGTACGTACGTACGT")])
    return g


class TestGetRegion:
    def test_basic_extraction(self, genome):
        assert genome.get_region("chr1", 1, 4) == "ACGT"

    def test_single_base(self, genome):
        assert genome.get_region("chr1", 1, 1) == "A"

    def test_full_sequence(self, genome):
        assert genome.get_region("chr1", 1, 20) == "ACGTACGTACGTACGTACGT"

    def test_reverse_complement(self, genome):
        # RC of ACGT is ACGT (palindrome)
        assert genome.get_region("chr1", 1, 4, strand="-") == "ACGT"
        # RC of ACG is CGT
        assert genome.get_region("chr1", 1, 3, strand="-") == "CGT"

    def test_out_of_bounds_raises(self, genome):
        with pytest.raises(ValueError):
            genome.get_region("chr1", 0, 5)
        with pytest.raises(ValueError):
            genome.get_region("chr1", 1, 21)

    def test_unknown_chromosome_returns_none(self, genome):
        assert genome.get_region("chrX", 1, 5) is None

    def test_start_gt_end_swaps(self, genome):
        # get_region normalizes start > end
        assert genome.get_region("chr1", 4, 1) == "ACGT"


class TestGenomeInfo:
    def test_has_chromosome(self, genome):
        assert genome.has_chromosome("chr1") is True
        assert genome.has_chromosome("chrX") is False

    def test_get_chromosome_size(self, genome):
        assert genome.get_chromosome_size("chr1") == 20
        assert genome.get_chromosome_size("chrX") is None

    def test_get_genome_size(self, genome):
        assert genome.get_genome_size() == 20

    def test_get_chromosome_names(self, genome):
        assert genome.get_chromosome_names() == ["chr1"]

    def test_get_base(self, genome):
        assert genome.get_base("chr1", 1) == "A"
        assert genome.get_base("chr1", 20) == "T"
        assert genome.get_base("chr1", 0) is None
        assert genome.get_base("chr1", 21) is None
        assert genome.get_base("chrX", 1) is None


class TestLoadFasta:
    def test_load_fasta(self, tmp_path):
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">contig1\nACGTACGT\n>contig2\nTTTTAAAA\n")
        g = Genome()
        names = g.load_fasta(str(fasta))
        assert names == ["contig1", "contig2"]
        assert g.get_region("contig1", 1, 4) == "ACGT"
        assert g.get_region("contig2", 5, 8) == "AAAA"

    def test_load_fasta_gz(self, tmp_path):
        import gzip
        fasta_gz = tmp_path / "test.fasta.gz"
        with gzip.open(str(fasta_gz), "wt") as f:
            f.write(">chr1\nGGGGCCCC\n")
        g = Genome()
        names = g.load_fasta(str(fasta_gz))
        assert "chr1" in names
        assert g.get_region("chr1", 1, 4) == "GGGG"


class TestIntMode:
    def test_int_mode_basic(self):
        g = Genome(storage="int")
        g.load_sequences([("chr1", "ACGTN")])
        # int mode reverse just reverses (no complement)
        assert g.get_region("chr1", 1, 5, strand="+") == "ACGTN"
        assert g.get_region("chr1", 1, 5, strand="-") == "NTGCA"
