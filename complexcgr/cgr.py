"From original work: CGR for gene structure"
from itertools import product
from tqdm import tqdm
from typing import Dict, Optional 
from collections import defaultdict, namedtuple
from PIL import Image
import numpy as np

# coordinates for x+iy
Coord = namedtuple("Coord", ["x","y"])

# coordinates for a CGR encoding
CGRCoords = namedtuple("CGRCoords", ["N","x","y"])

# coordinates for each nucleotide in the 2d-plane
DEFAULT_COORDS = dict(A=Coord(1,1),C=Coord(-1,1),G=Coord(-1,-1),T=Coord(1,-1))

class CGR: 
    "Chaos Game Representation for DNA"
    def __init__(self, coords: Optional[Dict[chr,tuple]]=None):
        self.nucleotide_coords = DEFAULT_COORDS if coords is None else coords
        self.cgr_coords = CGRCoords(0,0,0)
    
    def nucleotide_by_coords(self,x,y):
        "Get nucleotide by coordinates (x,y)"
        # filter nucleotide by coordinates
        filtered = dict(filter(lambda item: item[1] == Coord(x,y), self.nucleotide_coords.items()))

        return list(filtered.keys())[0]

    def forward(self, nucleotide: str): 
        "Compute next CGR coordinates"
        x = (self.cgr_coords.x + self.nucleotide_coords.get(nucleotide).x)/2
        y = (self.cgr_coords.y + self.nucleotide_coords.get(nucleotide).y)/2
        
        # update cgr_coords
        self.cgr_coords = CGRCoords(self.cgr_coords.N+1,x,y)

    def backward(self,):
        "Compute last CGR coordinates. Current nucleotide can be inferred from (x,y)"
        # get current nucleotide based on coordinates
        n_x,n_y = self.coords_current_nucleotide()
        nucleotide = self.nucleotide_by_coords(n_x,n_y)

        # update coordinates to the previous one
        x = 2*self.cgr_coords.x - n_x
        y = 2*self.cgr_coords.y - n_y
        
        # update cgr_coords
        self.cgr_coords = CGRCoords(self.cgr_coords.N-1,x,y)

        return nucleotide

    def coords_current_nucleotide(self,):
        x = 1 if self.cgr_coords.x>0 else -1
        y = 1 if self.cgr_coords.y>0 else -1
        return x,y

    def encode(self, sequence: str):
        "From DNA sequence to CGR"
        # reset starting position to (0,0,0)
        self.reset_coords()
        for nucleotide in sequence:
            self.forward(nucleotide)
        return self.cgr_coords
    
    def reset_coords(self,):
        self.cgr_coords = CGRCoords(0,0,0)

    def decode(self, N:int, x:int, y:int)->str: 
        "From CGR to DNA sequence"
        self.cgr_coords = CGRCoords(N,x,y)
        
        # decoded sequence
        sequence = []
        
        # Recover the entire genome
        while self.cgr_coords.N>0: 
            nucleotide = self.backward()
            sequence.append(nucleotide)
        return "".join(sequence[::-1])      