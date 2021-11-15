import random
from complexcgr import FCGR 

def test_fcgr():
    "frecuency matrix CGR"
    fcgr = FCGR(k=1)
    fcgr.count_kmers("ACGT")

    assert fcgr.freq_kmer.get("A") == 1
    assert fcgr.freq_kmer.get("C") == 1
    assert fcgr.freq_kmer.get("G") == 1
    assert fcgr.freq_kmer.get("T") == 1
    
    prob = 0.25
    fcgr.kmer_probabilities("ACGT")
    assert fcgr.probabilities.get("A") == prob
    assert fcgr.probabilities.get("C") == prob
    assert fcgr.probabilities.get("G") == prob
    assert fcgr.probabilities.get("T") == prob

def test_savefig():
    fcgr = FCGR(k=8)
    # Generate a random sequence without T's
    seq = "".join(random.choice("ACG") for _ in range(300_000))
    chaos = fcgr(seq) # an array with the probabilities of each k-mer
    fcgr.save(chaos, path="img/ACG.jpg")

def test_savefig_16bits():
    fcgr = FCGR(k=8, bits=16)
    # Generate a random sequence without T's
    seq = "".join(random.choice("ACG") for _ in range(300_000))
    chaos = fcgr(seq) # an array with the probabilities of each k-mer
    fcgr.save(chaos, path="img/ACG_16bits.jpg")




