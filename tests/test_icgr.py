from complexcgr import iCGR

def test_cgr():
    "Encoding for each nucleotide"
    cgr = iCGR()

    cgr.encode("A")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == 1
    assert cgr.cgr_coords.y == 1

    cgr.encode("C")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == -1
    assert cgr.cgr_coords.y == 1

    cgr.encode("G")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == -1
    assert cgr.cgr_coords.y == -1

    cgr.encode("T")
    assert cgr.cgr_coords.N == 1
    assert cgr.cgr_coords.x == 1
    assert cgr.cgr_coords.y == -1


def test_cgr_decode(): 
    "Decoding each nucleotide"
    cgr = iCGR()

    assert cgr.decode(N=1,x=1,y=1) == "A"
    assert cgr.decode(N=1,x=-1,y=1) == "C"
    assert cgr.decode(N=1,x=-1,y=-1) == "G"
    assert cgr.decode(N=1,x=1,y=-1) == "T"
    assert cgr.decode(N=4,x=3,y=-9) == "ACGT"