# Utils

def reverse_complement(sequence:str):
    """Returns the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
                 'a': 't', 'c': 'g', 'g': 'c', 't': 'a',
                 'N': 'N', 'n': 'n'}
    
    return ''.join(complement.get(base, base) for base in reversed(sequence))
