"From original work: CGR for gene structure"
from itertools import product
from tqdm import tqdm
from typing import Dict, Optional 
from collections import defaultdict, namedtuple

from matplotlib import pyplot as plt
from matplotlib import cm
import pylab

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

    def forward(self, nucleotide: str): 
        "Compute next CGR coordinates"
        x = (self.cgr_coords.x + DEFAULT_COORDS.get(nucleotide).x)/2
        y = (self.cgr_coords.y + DEFAULT_COORDS.get(nucleotide).y)/2
        
        # update cgr_coords
        self.cgr_coords = CGRCoords(self.cgr_coords.N+1,x,y)

    def backward(self,):
        "Compute last CGR coordinates. Current nucleotide can be inferred from (x,y)"
        # get current nucleotide based on coordinates
        n_x,n_y = self.coords_current_nucleotide()

        # update coordinates to the previous one
        x = 2*self.cgr_coords.x - n_x
        y = 2*self.cgr_coords.y - n_y
        
        # update cgr_coords
        self.cgr_coords = Coord(self.cgr_coords.N-1,x,y)

    def coords_current_nucleotide(self,):
        x = 1 if self.cgr_coords.x>0 else -1
        y = 1 if self.cgr_coords.y>0 else -1

    def encode(self, sequence: str):
        "From DNA sequence to CGR"
        # reset starting position to (0,0,0)
        self.reset_coords()
        for nucleotide in tqdm(sequence):
            self.forward(nucleotide)
    
    def reset_coords(self,):
        self.cgr_coords = CGRCoords(0,0,0)

    def decode(self, N:int, x:int, y:int)->str: 
        "From CGR to DNA sequence"
        self.cgr_coords = CGRCoords(N,x,y)
        # Recover the entire genome
        while self.cgr_coords.N>0: 
            self.backward()
        
class FCGR(CGR): 
    """Frequency matrix CGR
    an (4**kx4**k) 2D representation will be created for a sequence an 
    n-long sequence, k represents the k-mer
    """

    def __init__(self, k: int):
        super().__init__()
        self.k = k # k-mer representation
        self.freq_kmer = defaultdict(int)
        self.probabilities = defaultdict(float)
        self.kmers = product("ACGT", repeat=self.k)
        

    def count_kmers(self, sequence: str): 
        # representativity of kmers
        for j,_ in enumerate(sequence):
            subseq = sequence[j:j+self.k]
            if "N" not in subseq:
                self.freq_kmer[subseq] +=1
    
    def kmer_probabilities(self, sequence: str):
        N=len(sequence)
        for key, value in self.freq_kmer.items():
            self.probabilities[key] = float(value) / (N - self.k + 1)

    def __call__(self, sequence: str):
        "Frequency Matrix CGR"
        #https://towardsdatascience.com/chaos-game-representation-of-a-genetic-sequence-4681f1a67e14
        
        self.count_kmers(sequence)
        self.kmer_probabilities(sequence)

        array_size = int(2**self.k)
        fcgr = []
        for i in range(array_size):
            fcgr.append([0]*array_size)
    
        maxx = array_size
        maxy = array_size
        posx = 1
        posy = 1

        for kmer, prob in self.probabilities.items():
            for nb in kmer:
                if nb == "T":
                    posx += maxx / 2
                elif nb == "C":
                    posy += maxy / 2
                elif nb == "G":
                    posx += maxx / 2
                    posy += maxy / 2
                maxx = maxx / 2
                maxy /= 2
            
            fcgr[int(posy)-1][int(posx)-1] = prob
            maxx = array_size
            maxy = array_size
            posx = 1
            posy = 1
    
        return fcgr

    def plot(self, chaos):
        "Plot FCGR"
        pylab.title('Chaos game representation for {}-mers'.format(self.k))
        pylab.imshow(chaos, interpolation='nearest', cmap=cm.gray_r)
        pylab.show()