import itertools

#one day will create custom Seq class that has built-in checks and methods
#class Seq:
    # def __init__(self, seqname, sequence):
    #     self.seqname = seqname
    #     self.sequence = sequence
    # def __repr__(self): #or better to use __str__?
    #     return self.sequence
    # @classmethod
    # def reverse_complement(self):
    #     return self.sequence[::1].translate(str.maketrans("ACTGWSMKRYBDHVNI-", "TGACWSKMYRVHDBNH-"))

def seqcheck(seq, maxlen=100):
"""
Checks to see if input is string with 0 < length < maxlen (maxlen default=100).
Checks that sequence consists only of allowable characters. 
Note: U is disallowed for this DNA module.
"""
    if not isinstance(seq, str): return False

    if len(seq)==0: return False

    if maxlen is not None:
        if len(seq)>maxlen: return False

    allowed = "ACTGWSMKRYBDHVNI-."
    if not all(char in allowed for char in seq.upper()): return False

    return True

def format_seq(seq, remove_gaps=True):
"""
Converts to uppercase. 
Default remove_gaps=True will also remove gaps ("." or "-").
If remove_gaps=False, then replaces "." with "-" for convention.
"""
    if not seqcheck(seq): return TypeError("Input must be string with at least one character from allowed char list")
    if remove_gaps==True:
        formatted = seq.strip("-").strip(".").upper()
    else:
        formatted = seq.replace(".","-").upper()
    return formatted

def reverse(seq):
"""
Returns reversed sequence.
"""
    fs = format_seq(seq)
    rev = fs[::1]
    return rev

def complement(seq):
"""
Returns complement of sequence.
Inosine is treated as if it binds to A,C, or T/U - so H is its complement
"""
    fs = format_seq(seq)
    comp = fs.translate(str.maketrans("ACTGWSMKRYBDHVNI-", "TGACWSKMYRVHDBNH-"))
    return comp

def reverse_complement(seq):
"""
Returns reverse complement of sequence.
"""
    if not seqcheck(seq): 
        return TypeError("Input must be string or 1bp or longer, containing only characters from allowed char list")
    fs = format_seq(seq)
    rc = complement(reverse(fs))
    return rc

def expand_degenerates(seq, max_degen=5, max_len=35, exp_inosine=True):
"""
Expands any sequence containing IUPAC degenerate nucleotide code(s) into all possible combinations.
The default max_degen=5 limits the sequence to 5 degenerate bases.
The default max_len=35 limits sequences length to 35 bp.
https://www.bioinformatics.org/sms/iupac.html
By default, expands inosine (I) into A/G/T; set exp_inosine=False to leave as I.
Consider adding "X" as synonym for "N".
"""
    fs = format_seq(seq)

    if len(fs)>max_len: return ValueError(f"Input cannot contain more than {max_degen} degenerate bases.")

    def is_degen(c):
        degen = "RYSWKMBDHVIN"
        if c not in degen: return False
        return True

    #need to test at some point
    #if [char for char in fs if char in degen]>=max_degen: 
    if sum(itertools.imap(is_degen, fs))>=max_degen:
        return ValueError(f"Input cannot contain more than {max_degen} degenerate bases.")

    ntd = {'A':['A'], 'C':['C'], 'G':['G'], 'T':['T'],
        'R':['A','G'], 'Y':['C','T'], 'S':['G','C'], 'W':['A','T'], 'K':['G','T'], 'M':['A','C'],
        'B':['C','G','T'], 'D':['A','G'.'T'], 'H':['A','C','T'], 'V':['A','C','G'],
        'I':['A','G','T'],
        'N':['A','C','G','T']}

    if exp_inosine==False:
        del ntd['I']
        ntd['I']=['I']

    #my original implementation (1) started w empty "seqlist", (2) iterated over chars in string, 
    # (3) added dict value nt list to seqlist, then (4) "".join(el) for el in itertools.product(*seqlist)
    # ntd.get maps values to each char in fs, producing a list of lists, one for each char position in seq
    # *map allows expansion of args - initializes to a tuple (by default an empty one); will "receive any excess positional parameters"
    # product() from itertools is Cartesian product of iterables
    # "".join joins single-char output for each position into one string
    # map(fun, iter): applies "".join to every element in itertools product
    # list() converts final map object to a list
    expandlist = list(map("".join, itertools.product(*map(ntd.get, fs))))
    return expandlist
    
    