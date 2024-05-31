from fastahandle import fasta_check
from fastahandle import fasta2dict

#could probably just make this an extra argument flag for fasta_check
def validate_alignment(alignment, eqlength=True):
"""
Checks SIX conditions:
    1. File has valid fasta extension
    2. Filepath exists and can be opened
    3. All headers start with ">"
    4. No duplicate headers in fasta
    5. All seqs contain only valid characters
    ____________
    6. More than one sequence
    7. All seqs are same length (with default eqlength=True)
"""
    if not fasta_check(alignment): return False
    fdict = fasta2dict(alignment)

    if len(fdict)==1: return False

    if eqlength=True:
        lengthlist = [len(seq) for seq in fdict.values()]
        if set(lengthlist) != 1: return False

    return True




