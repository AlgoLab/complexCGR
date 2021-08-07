from . import complexCGR
import matplotlib.pyplot as plt
from itertools import product
from collections import defaultdict
from tqdm import tqdm
from PIL import Image
import numpy as np 


class complexFCGR(complexCGR): 
    """Circular density plot based on CGR"""

    def __init__(self, k: int):
        super().__init__()
        self.k = k # k-mer representation
        self.kmers = product("ACGT", repeat=self.k) # all kmers of length k
        self.freq_kmer = None # dict to save representativity of each kmer
        self.probabilities = None # dict to save probabilities/density for each kmer
        self.fig = None # to save matlotlib figure and then save as image
    
    def __call__(self, sequence: str, w=1):
        self.count_kmers(sequence)
        self.kmer_probabilities(sequence)
        self.plot(w)
    
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

    def plot(self, w: int = 1):        
        "Given a FCGR, plot it in grayscale"
        ax = plt.subplot(111, polar=True)
        center, bottom, width, height = self.compute_input_plot()
        
        # scale width by 'w'
        width = [w*_ for _ in width]
        print("generating plot")
        ax.bar(x=center, # center of the angle
                width=width, # width of the angle
                bottom=bottom, # lowest value
                height=height, # highest value 
            )

        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)

        return ax.figure

    def save(self, ccgr, path: str):
        "save complexFCGR as image"
        # get figure
        fig = self.plot() 

        # transform figure to Image
        img = self.fig2img(fig)
        
        # save image 
        img.save(path)

    @staticmethod
    def fig2img(fig):
        "Convert a Matplotlib figure to a PIL Image and return it"
        #https://stackoverflow.com/questions/57316491/how-to-convert-matplotlib-figure-to-pil-image-object-without-saving-image
        import io
        buf = io.BytesIO()
        fig.savefig(buf)
        buf.seek(0)
        img = Image.open(buf)
        return img

    
    def compute_input_plot(self,):
        "Compute input for plot CFCGR"
        delta = 2*np.pi/4**self.k # angle between consecutive roots

        center = []
        width  = []
        height = []
        bottom = []

        for kmer in tqdm(self.probabilities): 
            h = self.probabilities.get(kmer,0.0) # height (density) of the bar in the circle
            theta = 2*self.encode(kmer).k*np.pi/4**self.k # angle of k-esim root
            c = theta + delta/2 # center of the angle
            center.append(c)
            bottom.append(0) 
            width.append(50*delta)
            height.append(h)

        return center, bottom, width, height