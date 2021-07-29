from complexcgr import CGR

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

def test_cgr_decode(): 
    "Decoding each nucleotide"
    cgr = CGR()

    assert cgr.decode(N=1,x=0.5,y=0.5) == "A"
    assert cgr.decode(N=1,x=-0.5,y=0.5) == "C"
    assert cgr.decode(N=1,x=-0.5,y=-0.5) == "G"
    assert cgr.decode(N=1,x=0.5,y=-0.5) == "T"