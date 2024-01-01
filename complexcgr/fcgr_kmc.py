from . import FCGR

import gzip
import numpy as np

class FCGRKmc(FCGR):
    """
    Create FCGR with the option of using canonical kmers from KMC output
    """
    def __init__(self, k: int, use_canonical_kmers: bool=False):
        super().__init__(k, use_canonical_kmers)
        self.k = k # k-mer representation
        self.use_canonical_kmers = use_canonical_kmers

    def __call__(self, path_kmc_output):
        "Given a path to a kmc output file, return the FCGR using canonical kmers as an array"
        # Create an empty array to save the FCGR values
        array_size = int(2**self.k)
        fcgr = np.zeros((array_size,array_size))
        
        if str(path_kmc_output).endswith(".txt"):
            with open(path_kmc_output) as fp:
                for line in fp:
                    kmer, freq = line.split("\t")
                    pos_x, pos_y = self.kmer2pixel[kmer]
                    fcgr[int(pos_x)-1,int(pos_y)-1] += int(freq)
        else:
            with gzip.open(path_kmc_output,'rt') as f:
                for line in f:
                    line = line.strip()
                    kmer, freq = line.split()
                    pos_x, pos_y = self.kmer2pixel[kmer]
                    fcgr[int(pos_x)-1,int(pos_y)-1] += int(freq)
                    
        return fcgr