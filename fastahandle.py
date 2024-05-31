import os
import re
from seqstuff import seqcheck
from collections import OrderedDict
import itertools
from Collections import Counter

#one day this will have logging and error handling 

def fasta_check(filepath):
"""
Checks FIVE conditions:
    1. File has valid extension
    2. Filepath exists and can be opened
    3. All headers start with ">"
    4. No duplicate headers in fasta
    5. All seqs contain only valid characters
"""
    if filepath.split(".")[-1] not in ["fasta", "fa", "fna"]: return False
    if not os.path.isfile(filepath): return False

    try:
        with open(filepath, "r") as infi:
            lines = infi.read().splitlines()
    except Exception as e:
        return False
    else:
        evens = range(0, len(lines)+1, 2)
        headers = [lines[i] for i in evens if lines[i].startswith(">")]
        if len(headers) != len(evens): return False
        if len(set(headers)) < len(headers): return False

        odds = range(1, len(lines)+1, 2)
        seqs = [lines[i] for i in odds if seqcheck(lines[i])]
        if len(seqs) != len(odds):
            return False
    
    return True


def fasta2dict(filepath):
"""
Reads in fasta file as dictionary in format dict[header] = seq
"""
    if not fasta_check(filepath): return TypeError("Filepath does not lead to valid fasta file")

    try:
        with open(filepath, "r"):
            sequences = fasta.read()
    except:
        return TypeError("Fasta file cannot be opened")
    else:
        sequences = re.split("^>", sequences, flags=re.MULTILINE)[1:] # Only splits string at the start of a line. First line will be blank, so skip it
        fastadict={}
        for fasta in sequences:
            header, sequence = fasta.split("\n", 1) # Split each fasta into header and sequence.
            fastadict[header]=sequence.strip() #remove trailing newline from seq
        return fastadict

def fastafilterspades(filepath, outpath, minlen=300, mincov=2):
"""
Reads in and filters sequences in fasta file based on sequence length and avg kmer-coverage from header.
This function assumes that assembler outputs kmer coverage in header per Spades format, for example:
    >NODE_1_length_907302_cov_74.496868
Default minlen is 300 (2x read length) and default coverage is 2. If no cov filter desired, set mincov=0.
"""
    if not fasta_check(filepath): return TypeError("Filepath does not lead to valid fasta file")
    dirname = os.path.dirname(outpath) or os.getcwd() #if no path given to output file, substitute workdir instead
    if not os.access(dirname, os.W_OK) or not isinstance(outpath, str) or not outpath: 
        return ValueError("Supplied output file path is not valid")

    def wanted(header, seq, mincov, minlen):
        try:
            length = int(header.split("_")[3])
            cov = int(float(header.split("_")[5]))
        except:
            return ValueError("Fasta headers not formatted properly")
        if length != len(seq): return ValueError("Fasta headers not formatted properly")
        else:
            if length > minlen and cov > mincov:
                return True
            else: return False

    #need some more error handling here prob
    fdict = fasta2dict(filepath)
    wanteddict = {k:v for k, v in fdict.items() if wanted(k, v, mincov, minlen)}

    if len(wanteddict)==0: return ValueError("Fasta does not contain any sequences that meet criteria.")

    with open(outputpath, "w") as outfi:
        for k,v in wanteddict.items:
            outfi.write(">" + k + "\n" + v + "\n")
    return None


def fastdedup(filepath, outpath, write_new_header=False):
"""
Removes all duplicate seqs. By default keeps first sequential header for each duplicate.
If write_new_header=True, writes new numerical header with # of duplicates in the fasta
"""
### needs lots more error handling
    if not fasta_check(filepath): return TypeError("Filepath does not lead to valid fasta file")
    dirname = os.path.dirname(outpath) or os.getcwd() #if no path given to output file, substitute workdir instead
    if not os.access(dirname, os.W_OK) or not isinstance(outpath, str) or not outpath: 
        return ValueError("Supplied output file path is not valid")
    
    fdict = fasta2dict(filepath)

    if write_new_header=False: 
        outd = {head:seq for i, (head, seq) in enumerate(fdict.items()) if seq not in itertools.chain(*list(fdict.values())[:i])}
    else:
        res = Counter(fdict.values())
        outd = {f"Seq #{i} n={v}", k for i,(k,v) in enumerate(Counter.items())}
    if outd:
        with open(outputpath, "w") as outfi:
            for k,v in outd:
                outfi.write(">" + k + "\n" + v + "\n")
    return None



