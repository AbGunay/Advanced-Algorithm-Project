
# coding: utf-8

# In[100]:


# k mers motifs
# t - the numbers of dna strings
# dna - the array of strings
from collections import defaultdict
import operator
import time

def greedy_motif_search(dna, k):
    best_motifs = [d[0:3] for d in dna]
    first_str = dna[0]
    motif = []
    for i in range(len(first_str)-2):
        motif.append(first_str[i:i+k])
        for j in range(1, len(dna)):
            _profile = profile(motif,k)
            motif.append(probable_motif(_profile, dna[j],k))
        if score(motif) < score(best_motifs):
            best_motifs = motif
        motif = []
    return best_motifs

def profile(motif, k):
    profile = []
    l = len(motif)
    for i in range(k):
        s = {"A":0,"C":0,"G":0,"T":0}
        for j in range(len(motif)):
            if motif[j][i] == "A":
                s["A"] += 1.0
            elif motif[j][i] == "C":
                s["C"] += 1.0
            elif motif[j][i] == "G":
                s["G"] += 1.0
            else:
                s["T"] += 1.0
        for key, value in s.iteritems():
            s[key] = float("{0:.1f}".format(value/l))
        profile.append(s)    
    return profile
            

def probable_motif(profile, sequence, k):
    max = 0
    for s in range(len(sequence)-2):
        m = sequence[s:s+k]
        v = 1
        for p in profile:
            i = 0
            for key, value in p.iteritems():
                if m[i] == key:
                    v *= value  
                    break
            i += 1
        if v >= max:
            probable_motif = m
            max = v
    return probable_motif

def hamming(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def score(motifs):
    d = defaultdict(int)
    for motif in motifs:
        d[motif] += 1
    max_motif = max(d.iteritems(), key=operator.itemgetter(1))[0]
    s = 0
    for motif in motifs:
        s += hamming(max_motif, motif)
    return s

dna = ["GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG", "GGCGTTCAGGCA","AAGAATCAGTCA","CAAGGAGTTCGC","CACGTCAATCAC","CAATAATATTCG"]

start = time.time()
motif_search = greedy_motif_search(dna,3)
end = time.time()
print motif_search
print end - start


