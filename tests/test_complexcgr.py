from complexcgr import (
    __version__,
    complexCGR,
)

def test_version():
    assert __version__ == '0.7.0'

def test_complexCGR(): 
    ccgr = complexCGR()
    
    # nucleotide A
    encode = ccgr.encode("A")
    assert encode.k == 0
    assert encode.N == 1 


    # nucleotide C
    encode = ccgr.encode("C")
    assert encode.k == 1
    assert encode.N == 1


    # nucleotide G
    encode = ccgr.encode("G")
    assert encode.k == 2
    assert encode.N == 1


    # nucleotide T
    encode = ccgr.encode("T")
    assert encode.k == 3
    assert encode.N == 1

def test_complexCGR_decode():
    ccgr = complexCGR()

    # nucleotide A
    decode = ccgr.decode(k=0,N=1)
    assert decode == "A"


    # nucleotide C
    decode = ccgr.decode(k=1,N=1)
    assert decode == "C"


    # nucleotide G
    decode = ccgr.decode(k=2,N=1)
    assert decode == "G"


    # nucleotide T
    decode = ccgr.decode(k=3,N=1)
    assert decode == "T"

def test_encode_decode_compatibility():
    ccgr = complexCGR()

    seq = "ACGT"
    encode = ccgr.encode(seq)
    decode = ccgr.decode(encode.k, encode.N)
    assert seq == decode 