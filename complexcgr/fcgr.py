from . import CGR
from PIL import Image
from itertools import product
from collections import defaultdict
import numpy as np 

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
        len_seq = len(sequence)
        for j,_ in enumerate(sequence):
            if j+self.k <= len_seq:
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