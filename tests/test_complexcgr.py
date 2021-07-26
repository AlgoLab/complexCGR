from complexcgr import (
    __version__,
    CGR,
    FCGR,
)

def test_version():
    assert __version__ == '0.1.0'

def test_cgr():
    "Encoding for each nucleotide"
    cgr = CGR()

    cgr.encode("A")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == 0.5
    assert cgr.cgr_coords.y == 0.5

    cgr.encode("C")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == -0.5
    assert cgr.cgr_coords.y == 0.5

    cgr.encode("G")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == -0.5
    assert cgr.cgr_coords.y == -0.5

    cgr.encode("T")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == 0.5
    assert cgr.cgr_coords.y == -0.5

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