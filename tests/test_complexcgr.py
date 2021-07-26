from complexcgr import (
    __version__,
    CGR,
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