"""
GenomeHandler
Minimal, fast, and clean Python library for genomes and gene annotations.
"""

from .genome_db import GenomeDB
from .genome_class import Genome
from .gene_class import Gene

__all__ = ["GenomeDB", "Genome", "Gene"]

try:
    from importlib.metadata import version as _v
    __version__ = _v("genomehandler")
except Exception:
    __version__ = "0.0.0+local"