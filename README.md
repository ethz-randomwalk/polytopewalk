<p align="center">
  <img src="https://raw.githubusercontent.com/ethz-randomwalk/polytopewalk/main/docs/logo1.png" width="1000" object-fit = "cover">
</p>

![Python](https://img.shields.io/pypi/pyversions/polytopewalk.svg)
![PyPI](https://img.shields.io/pypi/v/polytopewalk)
![ciwheels](https://github.com/ethz-randomwalk/polytopewalk/actions/workflows/ciwheels.yml/badge.svg?branch=main)

# PolytopeWalk
**PolytopeWalk** is a `C++` library for running MCMC sampling algorithms to generate samples from a uniform distribution over a polytope with a `Python` interface. It handles preprocessing of the polytope (Facial Reduction algorithm) and initialization as well. Current implementations include the Dikin Walk, John Walk, Vaidya Walk, Ball Walk, Lee Sidford Walk, and Hit-and-Run in both the full-dimensional formulation and the sparse constrained formulation. For documentation on all functions/methods, please visit our webpage: https://polytopewalk.readthedocs.io/en/latest/ and read our paper on arXiv here: https://arxiv.org/abs/2412.06629. Finally, for example inputs and outputs, please visit the examples folder, which includes code to uniformly sample from both real-world polytopes from the `Netlib` dataset and structured polytopes.

## Code Structure

<p align="center">
  <img src="https://raw.githubusercontent.com/ethz-randomwalk/polytopewalk/main/docs/code_design1.png" width="750" object-fit = "cover">
</p>

## Implemented Algorithms
Let `d` be the dimension of the polytope, `n` be the number of boundaries, and `R/r` be where the convex body contains a ball of radius `r` and is mostly contained in a ball of radius `R`. We implement the following 6 MCMC sampling algorithms for uniform sampling over polytopes.  

| Name      | Mixing Time | Author |
| ------------ | ----------------- | ------------------- |
| `Ball Walk`   | $O(d^2R^2/r^2)$        | [Vempala (2005)](https://faculty.cc.gatech.edu/~vempala/papers/survey.pdf)       |
| `Hit and Run`   | $O(d^2R^2/r^2)$         | [Lovasz (1999)](https://link.springer.com/content/pdf/10.1007/s101070050099.pdf)         |
| `Dikin Walk`   | $O(nd)$         | [Sachdeva and Vishnoi (2015)](https://arxiv.org/pdf/1508.01977)     |
| `Vaidya Walk`   | $O(n^{1/2}d^{3/2})$        |   [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html)       |
| `John Walk`   | $O(d^{2.5})$        | [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html)           |
| `Lee Sidford Walk`   | $\tau(d^{2})$         | [Laddha et al. (2019)](https://arxiv.org/abs/1911.05656)  (conjectured, proof incomplete)         |

For each implemented algorithm, we provide the full-dimensional formulation and the sparse constrained formulation. Each polytope can be expressed from 1 formulation to the other. The main benefit of utilizing the constrained formulation is that it maintains sparse operations in A, ensuring scalability in higher dimensions. Many of the `netlib` dataset sparse polytopes are represented in this formulation. The formulations are specified below. 

In the full-dimensional formulation with dense matrix A ($n$ x $d$ matrix) and vector b ($n$ dimensional vector), we specify the following: 

```math
\mathcal{K}_1 = \{x \in \mathbb{R}^{d} | Ax \le b\}
```

where the polytope is specified with $n$ constraints.

In the constrained formulation with sparse matrix A ($n$ x $d$ matrix) and vector b ($n$ dimensional vector), we specify the following: 

```math
\mathcal{K}_2 = \{x \in \mathbb{R}^{d} | Ax = b, x \succeq_k 0\}
```

where the polytope is specified with $n$ equality constraints and $k$ coordinate-wise inequality constraints. 

In **PolytopeWalk**, we implement the MCMC algorithms in both the dense, full-dimensional and the sparse, constrained polytope formulation. 


## Installation

### Dependencies
**PolytopeWalk** requires:
- Python (>= 3.9)
- NumPy (>= 1.20)
- SciPy (>= 1.6.0)

### User installation
If you already have a working installation of NumPy and SciPy, the easiest way to install **PolytopeWalk** is using `pip`:
```bash
pip install -U polytopewalk
```


## Developer Installation Instructions 

### Important links
- Official source code repo: https://github.com/ethz-randomwalk/polytopewalk
- Download releases: https://pypi.org/project/polytopewalk/

### Install prerequisites
(listed in each of the operating systems)
- macOS: ``brew install eigen glpk``
- Linux:
    - Ubuntu ``sudo apt-get install -y libeigen3-dev libglpk-dev``
    - CentOS ``yum install -y epel-release eigen3-devel glpk-devel``
- Windows: ``choco install eigen -y``
    - Then, install winglpk from sourceforge

### Local install from source via pip
```bash
git clone https://github.com/ethz-randomwalk/polytopewalk.git
cd polytopewalk
pip install .
```


### Compile C++ from source (not necessary)
Only do this, if there is need to run and test C++ code directly. For normal users, we recommend only using the Python interface. 

Build with cmake
```bash
git clone https://github.com/ethz-randomwalk/polytopewalk.git && cd polytopewalk
cmake -B build -S . & cd build
make
sudo make install
```

## Examples
The `examples` folder provides examples of sampling from both sparse (constrained) and dense (full-dimensional) formulations of the MCMC sampling algorithms as well as testing convergence. We test our random walk algorithms on family of 3 structured polytopes and 3 polytopes from `netlib` for real-world analysis. The lines below show a quick demonstration of sampling from a polytope using a sparse MCMC algorithm. 
```python
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, csr_array
from polytopewalk.sparse import SparseDikinWalk

def generate_simplex(d):
    return np.array([1/d] * d), np.array([[1] * d]), np.array([1]), d, 'simplex'

x, A, b, k, name = generate_simplex(5)
sparse_dikin = SparseDikinWalk(r = 0.9)
dikin_res = sparse_dikin.generateCompleteWalk(10_000, x, A, b, k, burnin = 100, seed = 100)
```
We also demonstrate how to sample from a polytope in a dense, full-dimensional formulation. We additionally introduce the Facial Reduction algorithm, used to simplify the constrained polytope into the full-dimensional form. 
```python
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, csr_array
from polytopewalk.dense import DikinWalk, DenseCenter
from polytopewalk import FacialReduction

def generate_simplex(d):
    return np.array([1/d] * d), np.array([[1] * d]), np.array([1]), d, 'simplex'

fr = FacialReduction()
_, A, b, k, name = generate_simplex(5)
dikin = DikinWalk(r = 0.9)

polytope = fr.reduce(A, b, k, sparse = False)
dense_A = polytope.dense_A
dense_b = polytope.dense_b

dc = DenseCenter()
init = dc.getInitialPoint(dense_A, dense_b)

dikin_res = dikin.generateCompleteWalk(1_000, init, dense_A, dense_b, burnin = 100, seed = 100)
```

## Testing
The `tests` folder includes comprehensives tests of the Facial Reduction algorithm, Initialization, Weights from MCMC algorithms, and Sparse/Dense Random Walk algorithms in both Python and C++. Our Github package page comes with an automated test suite hooked up to continuous integration after push requests to the main branch. 

We provide instructions for locally testing **PolytopeWalk** in both Python and C++. For both, we must locally clone the repository (assuming we have installed the package already):

```bash
git clone https://github.com/ethz-randomwalk/polytopewalk.git
cd polytopewalk
```

### Python Testing
In addition to the requirements from the Developer Installation section, running this code requires a working version of Pandas. 

We can run the command:
```bash
python -m unittest discover -s tests/python -p "*.py"
```

### C++ Testing
As mentioned in the Developer Installation section, running this code requires a working version of Eigen and Glpk. 

First, we must compile the C++ code :
```bash
cmake -B build -S . && cd build 
make
```
Then, we can individually run the test files:
```bash
./tests/test_weights
./tests/test_fr
./tests/test_dense_walk
./tests/test_sparse_walk
./tests/test_init
```

## Community Guidelines

For those wishing to contribute to the software, please feel free to use the pull-request feature on our Github page, alongside a brief description of the improvements to the code. For those who have any issues with our software, please let us know in the issues section of our Github page. Finally, if you have any questions, feel free to contact the authors of this page at this email address: bys7@duke.edu.
