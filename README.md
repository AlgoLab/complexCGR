# complexCGR
`complexcgr` contains classes around the *Chaos Game Representation* for DNA sequences.

> version 0.5.0:  
A list of available classes and functionalities are listed below:
- [x] `CGR`  Chaos Game Representation: encodes a DNA sequence in 3 numbers $(N,x,y)$
  - [x] encode a sequence.
  - [x] recover a sequence from a CGR encoding.
- [x] `FCGR` Frequency Matrix CGR: representation as image for k-mer representativity, based on CGR.
  - [x] generate FCGR from an arbitrary n-long sequence.
  - [x] plot FCGR.
  - [x] save FCGR generated.
- [x] `iCGR` integer CGR: encodes a DNA sequence in 3 integers $(N,x,y)$. 
  - [x] encode a sequence
  - [x] recover a sequence from an iCGR encoding
- [x] `complexCGR`: encodes a DNA sequence in 2 integers $(k,N)$.
  - [x] encode a sequence
  - [x] recover a sequence from a complexCGR encoding
  - [x] plot sequence of complexCGR encodings 
- [x] `complexFCGR`: Frequency complexCGR: representation as image (circle) for k-mer representativity, based on complexCGR.
  - [x] generate complexFCGR from an arbitrary n-long sequence.
  - [x] plot complexFCGR.
  - [x] save complexFCGR generated.

## How to use
___
### 1. `CGR` Chaos Game Representation of DNA 
```python
from complexcgr import CGR

# Instantiate class CGR
cgr = CGR()

# encode a sequence
cgr.encode("ACGT")
# > CGRCoords(N=4, x=0.1875, y=-0.5625)

# recover a sequence from CGR coordinates
cgr.decode(N=4,x=0.1875,y=-0.5625)
# > "ACGT"
```

### 2. `FCGR` Frequency Matrix of Chaos Game Representation of DNA
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
| ![FCGR for a sequence without T's](img/CGA.jpg) |
|:--:|
|FCGR representation for a sequence without T's|


You can save the image with
```python
fcgr.save(chaos, path="img/ACG.jpg")
```

*Formats allowed are defined by PIL.*

```python
# Generate a random sequence without T's and lots of N's
seq = "".join(random.choice("ACGN") for _ in range(300_000))
chaos = fcgr(seq) # an array with the probabilities of each k-mer
fcgr.plot(chaos)
```


|![FCGR for a sequence without T's](img/CGAN.jpg)|
|:--:|
|FCGR representation for a sequence without T's and lot's of N's|


### 3. `iCGR` integer Chaos Game Representation of DNA 
```python
from complexcgr import iCGR

# Instantiate class CGR
icgr = iCGR()

# encode a sequence
icgr.encode("ACGT")
# > CGRCoords(N=4, x=3, y=-9)

# recover a sequence from CGR coordinates
icgr.decode(N=4,x=3,y=-9)
# > "ACGT"
```

### 4. `complexCGR` Complex Chaos Game Representation of DNA (complexCGR)

```python
from complexcgr import complexCGR

# Instantiate class CGR
ccgr = complexCGR()

# encode a sequence
ccgr.encode("ACGT")
# > CGRCoords(k=228,N=4)

# recover a sequence from complexCGR coordinates
ccgr.decode(k=228,N=4)
# > "ACGT"

```

### 5. `complexFCGR` Frequency Matrix of Complex Chaos Game Representation of DNA
Input for FCGR only accept sequences in $\{A,C,G,T,N\}$, but all $k$-mers that contains an $N$ 
will not be considered for the calculation of the frequency matrix CGR
```python
import random; random.seed(42)
from complexcgr import FCGR

# set the k-mer desired
cfcgr = complexFCGR(k=8) # 8-mers

# Generate a random sequence without T's
seq = "".join(random.choice("ACG") for _ in range(300_000))
fig = cfcgr(seq)

```
| ![FCGR for a sequence without T's](img/ACG-complexCGR.png) |
|:--:|
|complexFCGR representation for a sequence without T's|


You can save the image with
```python
cfcgr.save(fig, path="img/ACG-complexCGR.png")
```
*Currently the plot must be saved as png*


## Installation
___
```shell
pip install complexcgr
```

to update to the latest version
```shell
pip install complexcgr --upgrade
```