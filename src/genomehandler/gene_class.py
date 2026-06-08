
from dataclasses import dataclass, field

from .genome_utils import reverse_complement


@dataclass(slots=True)
class Gene:
    """Class representing a gene with its structure (exons, introns) and sequences."""

    gene_id: str
    gene_name: str
    chromosome: str
    start: int | None = None
    end: int | None = None
    strand: str | None = None
    transcripts: list[str] = field(default_factory=list)
    exons: dict[tuple[int, int], list[str]] = field(default_factory=dict)
    introns: dict[tuple[int, int], list[str]] = field(default_factory=dict)
    cdna_sequence: str | None = None
    protein_sequence: str | None = None
    protein_id: str | None = None
    product: str | None = None

    def add_exon(self, start: int, end: int, isoform: str | None = None):
        """Add an exon."""
        if isoform is None:
            isoform = self.gene_id

        if isoform not in self.transcripts:
            self.transcripts.append(isoform)

        if (start, end) not in self.exons:
            self.exons[(start, end)] = []

        if isoform not in self.exons[(start, end)]:
            self.exons[(start, end)].append(isoform)

    def calculate_introns(self) -> dict[tuple[int, int], list[str]]:
        """Calculate introns for each transcript based on exon positions."""
        transcript_exons: dict[str, list[tuple[int, int]]] = {}

        for (exon_start, exon_end), transcripts in self.exons.items():
            for transcript in transcripts:
                if transcript not in transcript_exons:
                    transcript_exons[transcript] = []
                transcript_exons[transcript].append((exon_start, exon_end))

        for transcript, exons in transcript_exons.items():
            sorted_exons = sorted(exons, key=lambda x: x[0])

            for i in range(len(sorted_exons) - 1):
                current_exon_end = sorted_exons[i][1]
                next_exon_start = sorted_exons[i + 1][0]

                if next_exon_start > current_exon_end + 1:
                    intron_start = current_exon_end + 1
                    intron_end = next_exon_start - 1

                    if (intron_start, intron_end) not in self.introns:
                        self.introns[(intron_start, intron_end)] = []

                    if transcript not in self.introns[(intron_start, intron_end)]:
                        self.introns[(intron_start, intron_end)].append(transcript)

        return self.introns

    def set_gene_start_end(self):
        """Set gene start and end coordinates based on exons."""
        if not self.exons:
            return None
        self.start = min(x[0] for x in self.exons.keys())
        self.end = max(x[1] for x in self.exons.keys())

    def extract_cdna_sequence(self, genome, transcript: str | None = None) -> str | None:
        """
        Extract coding sequence from genome for a specific transcript.
        If no transcript is specified, returns the sequence with all exons.

        Args:
            genome: Genome object with get_region method
            transcript: Optional transcript ID to extract specific isoform

        Returns:
            The cDNA sequence string or None if transcript not found
        """
        if not self.exons:
            return None

        relevant_exons: list[tuple[int, int]] = []
        if transcript:
            if transcript not in self.transcripts:
                return None

            for (start, end), isoforms in self.exons.items():
                if transcript in isoforms:
                    relevant_exons.append((start, end))

            if not relevant_exons:
                return None
        else:
            relevant_exons = list(self.exons.keys())

        # Sort exons by start position (ascending genomic order)
        relevant_exons.sort(key=lambda x: x[0])

        # Extract all exons on + strand, then RC the whole thing if minus
        cdna = ""
        for exon_start, exon_end in relevant_exons:
            cdna += genome.get_region(self.chromosome, exon_start, exon_end, strand="+")

        if self.strand == "-":
            cdna = reverse_complement(cdna)

        if not transcript:
            self.cdna_sequence = cdna

        return cdna

    # --- quick-access helpers ---
    def as_tuple(self, *, with_strand: bool = False) -> tuple:
        """Return (gene_id, gene_name, chromosome, start, end[, strand])."""
        if with_strand:
            return (self.gene_id, self.gene_name, self.chromosome, self.start, self.end, self.strand)
        return (self.gene_id, self.gene_name, self.chromosome, self.start, self.end)

    def region(self) -> str:
        """Return 'chrom:start-end' (or chrom:<unknown> if bounds missing)."""
        if self.start is None or self.end is None:
            return f"{self.chromosome}:<unknown>"
        return f"{self.chromosome}:{self.start}-{self.end}"

    def exons_as_tuples(self, *, sort: bool = True) -> list[tuple[int, int]]:
        """List exon intervals as (start, end)."""
        keys = list(self.exons.keys())
        if sort:
            keys.sort()
        return keys

    def introns_as_tuples(self, *, sort: bool = True) -> list[tuple[int, int]]:
        """List intron intervals as (start, end)."""
        keys = list(self.introns.keys())
        if sort:
            keys.sort()
        return keys

    def __call__(self) -> tuple:
        return self.as_tuple()

    def __str__(self) -> str:
        return f"Gene {self.gene_name} ({self.gene_id}): {self.chromosome}:{self.start}-{self.end} ({self.strand})"
