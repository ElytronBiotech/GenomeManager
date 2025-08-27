
from dataclasses import dataclass, field

@dataclass(slots=True)
class Gene:
    gene_id: str
    gene_name: str
    chromosome: str
    start: int | None = None
    end: int | None = None
    strand: str | None = None
    transcripts: list[str] = field(default_factory=list)
    exons: dict[tuple[int,int], list[str]] = field(default_factory=dict)
    introns: dict[tuple[int,int], list[str]] = field(default_factory=dict)
    cdna_sequence: str | None = None
    protein_sequence: str | None = None
    protein_id: str | None = None
    product: str | None = None

    """
    Class representing a gene with its structure (exons, introns) and sequences.
    """
    def __init__(self, gene_id, gene_name, chromosome, start=None, end=None, strand=None):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.chromosome = chromosome
        self.start = start  # 1-based
        self.end = end      # 1-based, inclusive
        self.strand = strand
        self.transcripts = []  # List of transcript IDs
        self.exons = {}     # Dict of (start, end) tuples for exons containing a list of isoforms
        self.introns = {}   # Dict of (start, end) tuples for introns containing a list of isoforms
        self.cdna_sequence = None
        self.protein_sequence = None
        self.protein_id = None
        self.product = None
    
    def add_exon(self, start:int, end:int, isoform:str = None):
        """Add an exon."""
        if isoform is None:
            isoform = self.gene_id
        
        # Add isoform to transcripts list if not already there
        if isoform not in self.transcripts:
            self.transcripts.append(isoform)
        
        if (start, end) not in self.exons:
            self.exons[(start, end)] = []
        
        if isoform not in self.exons[(start, end)]:
            self.exons[(start, end)].append(isoform)
    
    def calculate_introns(self):
        """Calculate introns for each transcript based on exon positions."""
        # Group exons by transcript
        transcript_exons = {}
        
        for (exon_start, exon_end), transcripts in self.exons.items():
            for transcript in transcripts:
                if transcript not in transcript_exons:
                    transcript_exons[transcript] = []
                transcript_exons[transcript].append((exon_start, exon_end))
        
        # Calculate introns for each transcript
        for transcript, exons in transcript_exons.items():
            # Sort exons by position
            sorted_exons = sorted(exons, key=lambda x: x[0])
            
            # Find introns between consecutive exons
            for i in range(len(sorted_exons) - 1):
                current_exon_end = sorted_exons[i][1]
                next_exon_start = sorted_exons[i+1][0]
                
                # Create intron if there's a gap
                if next_exon_start > current_exon_end + 1:
                    intron_start = current_exon_end + 1
                    intron_end = next_exon_start - 1
                    
                    # Add to introns dictionary
                    if (intron_start, intron_end) not in self.introns:
                        self.introns[(intron_start, intron_end)] = []
                    
                    if transcript not in self.introns[(intron_start, intron_end)]:
                        self.introns[(intron_start, intron_end)].append(transcript)
        
        return self.introns

    def set_gene_start_end(self):
        """Set gene start and end coordinates based on exons."""
        if not self.exons:
            return None
        #print(f'Changing gene start and end from {(self.start, self.end)} to:')
        self.start = min([x[0] for x in self.exons.keys()])
        self.end = max([x[1] for x in self.exons.keys()])
        #print(f'{(self.start, self.end)}')


    def extract_cdna_sequence(self, genome, transcript:str = None):
        """
        Extract coding sequence from genome for a specific transcript.
        If no transcript is specified, returns the sequence with all exons.
        
        Args:
            genome: Genome sequence string
            transcript: Optional transcript ID to extract specific isoform
            
        Returns:
            The cDNA sequence string or None if transcript not found
        """
        if not self.exons:
            return None
        
        # Filter exons by transcript if specified
        relevant_exons = []
        if transcript:
            # Check if the transcript exists
            if transcript not in self.transcripts:
                return None
                
            # Only include exons that contain this transcript
            for (start, end), isoforms in self.exons.items():
                if transcript in isoforms:
                    relevant_exons.append((start, end))
            
            # If no exons were found for this transcript, return None
            if not relevant_exons:
                return None
        else:
            # Include all exons
            relevant_exons = list(self.exons.keys())
        
        # Sort exons by start position
        relevant_exons.sort(key=lambda x: x[0])
        
        # Extract sequence
        cdna = ""
        for exon_start, exon_end in relevant_exons:
            # Adjust for 0-based indexing in Python
            cdna += genome.get_region(self.chromosome, exon_start, exon_end, strand=self.strand)  # <-- pass strand here
        
        # Only store in the instance variable if we're returning the full sequence
        if not transcript:
            self.cdna_sequence = cdna
        
        return cdna
    

    # --- quick-access helpers (most-used values) ---
    @property
    def product(self) -> str | None:
        # Prefer explicit product if you later choose to store it;
        # otherwise synthesize a reasonable label.
        if hasattr(self, "_product") and self._product:
            return self._product
        # fallback: prefer gene_name, else protein_id/name, else locus_tag/id
        return self.gene_name or getattr(self, "protein_id", None) or self.gene_id

    def as_tuple(self, *, with_strand: bool = False) -> tuple:
        """
        Return the common fields quickly.
        (gene_id, gene_name, chromosome, start, end[, strand])
        """
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

    # Optional: sugar so g() is the quick 4-tuple
    def __call__(self) -> tuple:
        return self.as_tuple()
        
    def __str__(self):
        return f"Gene {self.gene_name} ({self.gene_id}): {self.chromosome}:{self.start}-{self.end} ({self.strand})"