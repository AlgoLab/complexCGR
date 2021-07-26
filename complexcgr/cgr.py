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
        for nucleotide in sequence:
            self.forward(nucleotide)
        return self.cgr_coords
    
    def reset_coords(self,):
        self.cgr_coords = CGRCoords(0,0,0)

    def decode(self, N:int, x:int, y:int)->str: 
        "From CGR to DNA sequence"
        self.cgr_coords = CGRCoords(N,x,y)
        # Recover the entire genome
        while self.cgr_coords.N>0: 
            self.backward()
        return self.cgr_coords        

class FCGR(CGR): 
    """Frequency matrix CGR
    an (2**k x 2**k) 2D representation will be created for a 
    n-long sequence. 
    - k represents the k-mer.
    - 2**k x 2**k = 4**k the total number of k-mers (sequences of length k)
    """

    def __init__(self, k: int):
        super().__init__()
        self.k = k # k-mer representation
        self.kmers = product("ACGT", repeat=self.k)

    def count_kmers(self, sequence: str): 
        self.freq_kmer = defaultdict(int)
        # representativity of kmers
        for j,_ in enumerate(sequence):
            subseq = sequence[j:j+self.k]
            if "N" not in subseq:
                self.freq_kmer[subseq] +=1
    
    def kmer_probabilities(self, sequence: str):
        self.probabilities = defaultdict(float)
        N=len(sequence)
        for key, value in self.freq_kmer.items():
            self.probabilities[key] = float(value) / (N - self.k + 1)

    def pixel_position(self, kmer: str):
        "Get pixel position in the FCGR matrix for a k-mer"

        coords = self.encode(kmer)
        N,x,y = coords.N, coords.x, coords.y
        
        # Coordinates from [-1,1]² to [1,2**k]²
        np_coords = np.array([(x + 1)/2, (y + 1)/2]) # move coordinates from [-1,1]² to [0,1]²
        np_coords *= 2**self.k # rescale coordinates from [0,1]² to [0,2**k]²
        x,y = np.ceil(np_coords) # round to upper integer 

        # Turn coordinates (cx,cy) into pixel (px,py) position 
        # px = 2**k-cy+1, py = cx
        return 2**self.k-int(y)+1, int(x)

    def __call__(self, sequence: str):
        "Given a DNA sequence, returns an array with his FCGR"
        self.count_kmers(sequence)
        self.kmer_probabilities(sequence)
        # Create an empty array to save the FCGR values
        array_size = int(2**self.k)
        fcgr = np.zeros((array_size,array_size))

        # Assign probability to each box in the Frequency CGR
        for kmer, prob in self.probabilities.items():
            pos_x, pos_y = self.pixel_position(kmer)          
            fcgr[int(pos_x)-1,int(pos_y)-1] = prob
        return fcgr

    def plot(self, fcgr):
        "Given a FCGR, plot it in grayscale"
        img_pil = self.array2img(fcgr)
        return img_pil

    def save(self, fcgr, path: str):
        "Save image in grayscale for the FCGR provided as input"
        img_pil = self.array2img(fcgr)
        img_pil.save(path)
    
    @staticmethod
    def array2img(array):
        "Array to PIL image"
        m, M = array.min(), array.max()
        # rescale to [0,1]
        img_rescaled = (array - m) / (M-m) 
        
        # invert colors black->white
        img_array = np.ceil(255 - img_rescaled*255)
        img_array = np.array(img_array, dtype="uint8")
        
        # convert to Image 
        img_pil = Image.fromarray(img_array,'L')
        return img_pil