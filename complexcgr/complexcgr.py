from typing import Dict
from collections import namedtuple

# coordinates for x+iy
Coord = namedtuple("Coord", ["x","y"])

# coordinates for a CGR encoding
CGRCoords = namedtuple("CGRCoord", ["k","N"])

# coordinates for each nucleotide in the complex plane
DEFAULT_COORDS = dict(A=Coord(1,0),C=Coord(0,1),G=Coord(-1,0),D=Coord(0,-1))

class complexCGR: 
    "Complex Chaos Game Representation of DNA"
    def __init__(self, coords: Dict[chr,tuple]=DEFAULT_COORDS):
        self.nucleotide_coords = coords

    def forward(self,): 
        "compute next complexCGR coordinates"
        pass

    def backward(self,):
        "compute last complexCGR coordinates"
        pass
    