from complexcgr import (
    __version__,
    complexCGR,
)

def test_version():
    assert __version__ == '0.2.0'

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