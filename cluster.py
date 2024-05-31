#hammingclust, findcentroid
import sklearn.cluster as cluster
from seqstuff import seqcheck
from fastahandle import fasta2dict
from alignment import validate_alignment

def hammingdist(s1, s2):
"""
Finds pairwise Hamming distance between two seqs.
"""
    if not seqcheck(s1) and seqcheck(s2): return ValueError("Sequences not in valid format")

    return sum(el1 != el2 for el1 el2 in zip(s1, s2))

#matrix, numclust, matrix, gapopen, gapextend):
def kmeansclust(alignment, numclust, outdir):
"""
Uses Kmeans clustering to cluster an all-vs-all Hamming distance matrix.
Can also include option to use aligner scoring instead of Hamming dist (ex biopython pairwise2) 
outdir is the directory where all otuput cluster fastas will be written to.
"""
    if not validate_alignment(alignment): return ValueError("Alignment not in valid format")

    #validate outdir

    #get numbered list of header-seqs
    alndict = fast2dict(alignment)
    recs={}
    [recs[i]=k:v for i,(k,v) in enumerate(alndict.items())]

    #make empty distance matrix 
    #distm = [[0 for i in range(len(recs))] for j in range(len(recs))]

    #make dataframe with key names as index
    df = pd.DataFrame(columns=fastadict.keys(),index=fastadict.keys())
    #iterate through df and put distance metrics into cells
    for i,r in df.iterrows():
        for col in df.columns:
            query=fastadict[r.name]
            seqlength=len(query)
            ref=fastadict[col]
            hamming = hamming_distance(query,ref)
            hamdist=(seqlength-hamming)/seqlength #proportion of same-character positions to total sequence length
            df.at[r.name,col] = hamdist
        # ^ df.at????

    #double iterate over numbered list. 
    #for i in range(0, len(recs)):
    #    for j in range(0, len(recs)):
    #        seq1 = recs[i]
    #        seq2 = recs[j]
    #        hd = hammingdist(seq1, seq2)
    #        distm[i][j] = hd
    #df to matrix
    mat = df.values

    kmeans = sklearn.cluster.KMeans(numclust)
    results = kmeans.fit(distm)
    labels = results.labels_

    for i in len(df.index):
        header = ">" + labels[i] + " " + df.index[i] + "\n"
        seq = #enumerated sequence

    clusters = [[] for i in range(n_clusters)]

    for i in range(0, len(recs)):
        clusters[labels[i]].append(recs[i])
    for i in range(0, len(clusters)):

        
        


