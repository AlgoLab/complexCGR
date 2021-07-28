# complexCGR
`complexcgr` contains classes around the Game Representation for DNA sequences.
The available classes are: 
- [x] `CGR`  Chaos Game Representation: encodes a DNA sequence in 3 numbers $(N,x,y)$
- [x] `FCGR` Frequency Matrix CGR: representation as image for k-mer representativity based on CGR.
- [ ] `iCGR` integer CGR: encodes a DNA sequence in 3 integers $(N,x,y)$. 
- [x] `complexCGR`: encodes a DNA sequence in 2 integers $(k,N)$.

## How to use
___
### 1. Chaos Game Representation of DNA (CGR)
```python
from complexcgr import CGR

# Instantiate class CGR
cgr = CGR()
cgr_encode = cgr.encode("ACGT")
# > CGRCoords(N=4, x=0.1875, y=-0.5625)
```

### 2. Frequency Matrix of Chaos Game Representation of DNA (FCGR)
Input for FCGR only accept sequences in $\{A,C,G,T,N\}$, but all $k$-mers that contains an $N$ 
will not be considered for the calculation of the frequency matrix CGR
```python
import random; random.seed(42)
from complexcgr import FCGR

# set the k-mer desired -> (2**k,2**k)$ FCGR
fcgr = FCGR(k=8) # (256x256) image

# Generate a random sequence without T's
seq = "".join(random.choice("ACG") for _ in range(300_000))
chaos = fcgr(seq) # an array with the probabilities of each k-mer
fcgr.plot(chaos)
```
| ![CGR for a sequence without T's](img/CGA.jpg) |
|:--:|
|sequence without T's|

```python
# Generate a random sequence without T's and lots of N's
seq = "".join(random.choice("ACGN") for _ in range(300_000))
chaos = fcgr(seq) # an array with the probabilities of each k-mer
fcgr.plot(chaos)
```


|![CGR for a sequence without T's](img/CGAN.jpg)|
|:--:|
|sequence without T's and lot's of N's|

### 3. Complex Chaos Game Representation of DNA (complexCGR)

```python
# Instantiate class CGR
ccgr = complexCGR()
ccgr_encode = ccgr.encode("ACGT")
# > CGRCoords(k=228,N=4)
```


## Installation
___
```shell
pip install complexcgr
```

to update to the latest version
```shell
pip install complexcgr --upgrade
```