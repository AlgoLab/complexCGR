import numpy as np
from Bio import SeqIO
from collections import defaultdict
from pathlib import Path
from tqdm import tqdm
from typing import List, Union
from complexcgr import FCGR 

# input for FCGR Samples
_path_fastq = Union[str, Path] # path can be a string or a Path instance
_fastq = Union[_path_fastq,List[_path_fastq]] # either a single path or a list of paths to fastq files

class FCGRSamples(FCGR):

    def __init__(self, k: int, bits: int = 8):
        super().__init__(k, bits)

    def __call__(self, path_fastq: _fastq, consider_quality: bool = False):
        "Given a (list) of fastq files, return the FCGR matrix" 
        
        ## ---- New for reads ---- 
    
        # for each call, we start to count the frequencies
        self.freq_kmer = defaultdict(int)
        if consider_quality is True:
            self.qual_kmer = defaultdict(int)
        
        # transform to list to iterate
        path_fastq = path_fastq if type(path_fastq) is list else [path_fastq]
        
        # For each file, count kmers on their reads
        for path in path_fastq: 
            fastq = self.load_fastq(path)
            n_reads = self.count_reads(path)

            # count kmers in each read
            for sample in tqdm(fastq,total=n_reads,desc=f"Counting kmers on {str(Path(path).stem)}"):
                read = sample.seq
                if consider_quality is True:
                    qual = sample.letter_annotations["phred_quality"]
                    self.count_kmers_qualities(str(read),qual)
                else:
                    self.count_kmers(str(read))

        ## ---- End new for reads ----

        # Create an empty array to save the FCGR values
        if consider_quality is False:
            array_size = int(2**self.k)
            fcgr = np.zeros((array_size,array_size))
        else:
            array_size = int(2**self.k)
            fcgr = np.zeros((array_size,array_size,2))

        # Assign frequency to each box in the matrix
        for kmer, freq in self.freq_kmer.items():        
            pos_x, pos_y = self.kmer2pixel[kmer]
            
            if consider_quality is True:
                fcgr[int(pos_x)-1,int(pos_y)-1,0] = freq
                fcgr[int(pos_x)-1,int(pos_y)-1,1] = self.qual_kmer[kmer]
            else:
                fcgr[int(pos_x)-1,int(pos_y)-1] = freq

        return fcgr  if consider_quality is False else self.rescale_fcgr_qualities(fcgr)
        

    def load_fastq(self, path): 
        "Load a fastq file"
        fastq = SeqIO.parse(str(path), "fastq")
        return fastq

    def count_kmers(self, read: str): 
        # representativity of kmers
        last_j = len(read) - self.k + 1   
        kmers  = (read[i:(i+self.k)] for i in range(last_j))
        # count kmers in a dictionary
        list(self.count_kmer(kmer) for kmer in kmers)        

    @staticmethod
    def count_reads(path):
        "Count reads in a file to define progress bar"
        n_reads = 0
        with open(str(path)) as fp: 
            for line in fp: 
                if line.startswith("@"): 
                    n_reads += 1
        return n_reads

    def count_kmer(self, kmer):
        if "N" not in kmer:
            self.freq_kmer[kmer] += 1

    
    def count_kmer_quality(self, kmer, qmer):
        """Count kmers and qualities
        For quality of a kmer, the mean of qualities for each nucleotide in the kmer will be saved"""
        if "N" not in kmer:
            self.freq_kmer[kmer] += 1
            self.qual_kmer[kmer] += qmer.mean()

    def count_kmers_qualities(self,read,qual):
        "Count kmers and qualities"

        # representativity of kmers
        last_j = len(read) - self.k + 1   
        kmers  = (read[i:(i+self.k)] for i in range(last_j)) 
        qmers  = (np.array(qual[i:(i+self.k)]) for i in range(last_j)) # qualities for each kmer
        
        # count kmers and qualities in a dictionary
        list(self.count_kmer_quality(kmer,qmer) for (kmer,qmer) in zip(kmers,qmers))

    @staticmethod
    def rescale_fcgr_qualities(fcgr):
        "Divide each cumulate quality by the freq of its kmer"
        freqs, quals = fcgr[:,:,0], fcgr[:,:,1]
        fcgr[:,:,1] = np.divide(quals, freqs, out = np.zeros_like(quals), where = quals!=0)
        return fcgr