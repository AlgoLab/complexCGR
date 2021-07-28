from complexcgr.cgr import CGR
from typing import Dict, List
from collections import namedtuple

# coordinates for x+iy
Coord = namedtuple("Coord", ["x","y"])

# coordinates for a CGR encoding
CGRCoords = namedtuple("CGRCoord", ["k","N"])

class complexCGR:
    "Complex Chaos Game Representation of DNA"
    def __init__(self, nucleotide_order: List[str] = ["A","C","G","T"]):
        self.nucleotide_order = nucleotide_order
        self.cgr_coords = CGRCoords(0,0) # complexCGR coordinates 
    
    def id(self, nucleotide): 
        return self.nucleotide_order.index(nucleotide)
            
    def forward(self,nucleotide: str): 
        "compute next complexCGR coordinates"
        k = self.id(nucleotide)*4**(self.cgr_coords.N) + self.cgr_coords.k

        # update cgr_coords
        self.cgr_coords = CGRCoords(k,self.cgr_coords.N+1)

    def backward(self,):
        "compute last complexCGR coordinates"
        nucleotide = self.current_nucleotide()    
        k = self.cgr_coords.k - self.id(nucleotide)*4**(self.cgr_coords.N)

        # update cgr_coords
        self.cgr_coords = CGRCoords(k,self.cgr_coords.N-1)
    
    def encode(self, sequence: str): 
        "From DNA to complexCGR"
        self.reset_coords()
        for nucleotide in sequence: 
            self.forward(nucleotide)
        return self.cgr_coords

    def decode(self, k: int, N: int): 
        "From complexCGR to DNA"
        self.cgr_coords = CGRCoords(k,N)
        # Recover the entire genome
        while self.cgr_coords.N>0: 
            self.backward()
        return self.cgr_coords        

    def current_nucleotide(self,): 
        "Get current nucleotide based on k and N"
        alpha = self.k/self.N
        if alpha <0.25: 
            return self.nucleotide_order[0]
        elif alpha <0.5:
            return self.nucleotide_order[1]
        elif alpha < 0.75:
            return self.nucleotide_order[2]
        else:
            return self.nucleotide_order[3]

    def reset_coords(self,):
        self.cgr_coords = CGRCoords(0,0)
