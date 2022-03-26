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

    def __init__(self, k: int, bits: int = 8):
        super().__init__()
        self.k = k # k-mer representation
        self.kmers = list("".join(kmer) for kmer in product("ACGT", repeat=self.k))
        self.kmer2pixel = self.kmer2pixel_position()
        self.freq_kmer = defaultdict(int)
        self.bits = bits
        self.max_color = 2**bits-1

    def __call__(self, sequence: str):
        "Given a DNA sequence, returns an array with his FCGR"
        self.count_kmers(sequence)
       
        # Create an empty array to save the FCGR values
        array_size = int(2**self.k)
        fcgr = np.zeros((array_size,array_size))

        # Assign frequency to each box in the matrix
        for kmer, freq in self.freq_kmer.items():        
            pos_x, pos_y = self.kmer2pixel[kmer]
            fcgr[int(pos_x)-1,int(pos_y)-1] = freq
        return fcgr

    def count_kmer(self, kmer):
        if "N" not in kmer:
            self.freq_kmer[kmer] += 1

    def count_kmers(self, sequence: str): 
        self.freq_kmer = defaultdict(int)
        # representativity of kmers
        last_j = len(sequence) - self.k + 1   
        kmers  = (sequence[i:(i+self.k)] for i in range(last_j))
        # count kmers in a dictionary
        list(self.count_kmer(kmer) for kmer in kmers)

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

    def kmer2pixel_position(self,):
        kmer2pixel = dict()
        for kmer in self.kmers:
            kmer2pixel[kmer] = self.pixel_position(kmer)
        return kmer2pixel

    def plot(self, fcgr):
        "Given a FCGR, plot it in grayscale"
        img_pil = self.array2img(fcgr)
        return img_pil

    def save_img(self, fcgr, path: str):
        "Save image in grayscale for the FCGR provided as input"
        img_pil = self.array2img(fcgr)
        img_pil.save(path)
    
    def array2img(self, array):
        "Array to PIL image"
        m, M = array.min(), array.max()
        # rescale to [0,1]
        img_rescaled = (array - m) / (M-m) 
        
        # invert colors black->white
        img_array = np.ceil(self.max_color - img_rescaled*self.max_color)
        dtype = eval(f"np.int{self.bits}")
        img_array = np.array(img_array, dtype=dtype)
        
        # convert to Image 
        img_pil = Image.fromarray(img_array,'L')
        return img_pil