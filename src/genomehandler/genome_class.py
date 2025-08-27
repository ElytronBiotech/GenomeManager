import gzip
import numpy as np
from Bio import SeqIO

_ASCII_RC = np.full(256, ord('N'), dtype=np.uint8)  # default complement for unknowns
_ASCII_RC[ord('A')] = ord('T'); _ASCII_RC[ord('a')] = ord('t')
_ASCII_RC[ord('C')] = ord('G'); _ASCII_RC[ord('c')] = ord('g')
_ASCII_RC[ord('G')] = ord('C'); _ASCII_RC[ord('g')] = ord('c')
_ASCII_RC[ord('T')] = ord('A'); _ASCII_RC[ord('t')] = ord('a')

class Genome:
    """
    Memory-efficient genome representation.
    By default, stores sequences as ASCII bytes (np.uint8), 1 byte/base.
    """
    def __init__(self, *, storage: str = "bytes"):
        """
        storage:
          - "bytes": ASCII bytes (np.uint8), 1 byte/base (default)
          - "int":   small integer codes (A=1,C=2,G=3,T=4,N=5), also 1 byte/base
        """
        assert storage in {"bytes", "int"}
        self._storage = storage
        self.chromosomes: dict[str, np.ndarray] = {}

    # --- helpers for encoding/decoding ---
    def _encode(self, s: str) -> np.ndarray:
        if self._storage == "bytes":
            return np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        # "int" mode
        # vectorized mapping via lookup table
        lut = np.zeros(256, dtype=np.uint8)
        lut[ord('A')] = 1; lut[ord('C')] = 2; lut[ord('G')] = 3; lut[ord('T')] = 4; lut[ord('N')] = 5
        lut[ord('a')] = 1; lut[ord('c')] = 2; lut[ord('g')] = 3; lut[ord('t')] = 4; lut[ord('n')] = 5
        arr = np.frombuffer(s.encode("ascii"), dtype=np.uint8)
        return lut[arr]

    def _decode(self, arr: np.ndarray) -> str:
        if self._storage == "bytes":
            return arr.tobytes().decode("ascii")
        # "int" mode decode via lookup
        lut = np.array([ord('?'), ord('A'), ord('C'), ord('G'), ord('T'), ord('N')], dtype=np.uint8)
        ascii_arr = lut[arr]
        return ascii_arr.tobytes().decode("ascii")

    # --- public API ---
    def load_fasta(self, fasta_file: str):
        """
        Load genome from a FASTA file, storing each chromosome as a compact np.uint8 array.
        """
        self.chromosomes.clear()
        open_fn = gzip.open if fasta_file.endswith(".gz") else open
        with open_fn(fasta_file, "rt") as handle:
            for rec in SeqIO.parse(handle, "fasta"):
                self.chromosomes[rec.id] = self._encode(str(rec.seq))
        return list(self.chromosomes.keys())
    
    # Added this function to handle gbk integrated genomes
    def load_sequences(self, seqs: list[tuple[str, str]]):
        self.chromosomes.clear()
        for name, s in seqs:
            self.chromosomes[name] = self._encode(s)
        return list(self.chromosomes.keys())

    def has_chromosome(self, chrom_name: str) -> bool:
        return chrom_name in self.chromosomes

    def get_chromosome(self, chrom_name: str) -> np.ndarray | None:
        return self.chromosomes.get(chrom_name)

    def get_base(self, chrom_name: str, position: int) -> str | None:
        """
        1-based position access; returns a single-character string.
        """
        arr = self.chromosomes.get(chrom_name)
        if arr is None:
            return None
        i = position - 1
        if i < 0 or i >= arr.shape[0]:
            return None
        if self._storage == "bytes":
            return chr(int(arr[i]))
        # int mode
        return self._decode(arr[i:i+1])

    def get_region(self, chrom_name: str, start: int, end: int, *, strand: str = "+") -> str | None:
        """
        Return the sequence string for 1-based inclusive region [start, end].
        If strand == '-', reverse-complement is returned (case preserved for ASCII; int mode upper-case).
        """
        arr = self.chromosomes.get(chrom_name)
        if arr is None:
            return None
        # normalize
        if start > end:
            start, end = end, start
        i0 = start - 1
        i1 = end     # slicing end is exclusive
        if i0 < 0 or i1 > arr.shape[0]:
            raise ValueError(f"Region out of bounds: {chrom_name}:{start}-{end}")

        view = arr[i0:i1]  # zero-copy slice
        if strand == "+":
            return self._decode(view)

        # reverse-complement for '-' (only precise in bytes mode; in int mode we just reverse)
        if self._storage == "bytes":
            rc = _ASCII_RC[view][::-1]
            return rc.tobytes().decode("ascii")
        else:
            # int mode: reverse + decode (no complementing to keep mapping simple)
            return self._decode(view[::-1])
    
    def get_chromosome_size(self, chrom_name):
        """Get the length of a chromosome."""
        if chrom_name not in self.chromosomes:
            return None
        return len(self.chromosomes[chrom_name])
    
    def get_genome_size(self):
        """Get the total size of the genome."""
        return sum(len(seq) for seq in self.chromosomes.values())
    
    def get_chromosome_names(self):
        """Get a list of all chromosome names."""
        return list(self.chromosomes.keys())
    
    def __str__(self):
        """String representation of the genome info."""
        chrom_info = [f"{name} ({len(seq)} bp)\n" for name, seq in self.chromosomes.items()]
        chrom_info_str = "    ".join(chrom_info)
        return f"Genome with {len(self.chromosomes)} chromosomes:\n    {chrom_info_str}"