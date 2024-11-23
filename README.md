<p align="center">
  <img src="https://raw.githubusercontent.com/ethz-randomwalk/polytopewalk/main/docs/logo1.png" width="1000" object-fit = "cover">
</p>

![Python](https://img.shields.io/pypi/pyversions/polytopewalk.svg)
![PyPI](https://img.shields.io/pypi/v/polytopwalk)
![ciwheels](https://github.com/ethz-randomwalk/polytopewalk/actions/workflows/ciwheels.yml/badge.svg?branch=main)

# PolytopeWalk
**PolytopeWalk** is a `C++` library for running MCMC sampling algorithms to generate samples from a uniform distribution over a polytope with a `Python` interface. It handles preprocessing of the polytope and initialization as well. Current implementations include the Dikin Walk, John Walk, Vaidya Walk, Ball Walk, Lee Sidford Walk, and Hit-and-Run in both the full-dimensional formulation and the sparse constrained formulation. Code includes facial reduction and initialization algorithms for pre-processing as well. Sample code that samples from both real polytopes from a data set and artificial polytopes are shown in the Examples folder.

## Code Structure

<p align="center">
  <img src="https://raw.githubusercontent.com/ethz-randomwalk/polytopewalk/main/docs/code_design1.png" width="750" object-fit = "cover">
</p>

## Implemented Algorithms
Let `d` be the dimension of the polytope, `n` be the number of boundaries, and `R/r` be where the convex body contains a ball of radius `r` and is mostly contained in a ball of radius `R`. We implement the following 6 MCMC sampling algorithms for uniform sampling over polytopes.  

| Name      | Mixing Time | Author |
| ------------ | ----------------- | ------------------- |
| `Ball Walk`   | $\tau(d^2R^2/r^2)$        | [Vempala (2005)](https://faculty.cc.gatech.edu/~vempala/papers/survey.pdf)       |
| `Hit and Run`   | $\tau(d^2R^2/r^2)$         | [Lovasz (1999)](https://link.springer.com/content/pdf/10.1007/s101070050099.pdf)         |
| `Dikin Walk`   | $\tau(nd)$         | [Sachdeva and Vishnoi (2015)](https://arxiv.org/pdf/1508.01977)     |
| `Vaidya Walk`   | $\tau(n^{1/2}d^{3/2})$        |   [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html)       |
| `John Walk`   | $\tau(d^{2.5})$        | [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html)           |
| `Lee Sidford Walk`   | $\tau(d^{2})$         | [Laddha et al. (2019)](https://arxiv.org/abs/1911.05656)          |

## Installation

### Dependencies
polytopewalk requires:
- Python (>= 3.9)
- NumPy (>= 1.20)
- SciPy (>= 1.6.0)

### User installation
If you already have a working installation of NumPy and SciPy, the easiest way to install polytopewalk is using `pip`:
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
cmake -B build -S .
make
sudo make install
```

## Examples
The `examples` folder provides examples of sampling from both sparse and dense formulations of the MCMC sampling algorithms. We test our random walk algorithms on family of 3 structured polytopes and 3 real life polytopes from `netlib`. The lines below show a quick demonstration of sampling from a polytope using a sparse MCMC algorithm. 
```python
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, csr_array
from polytopewalk.sparse import SparseDikinWalk

def generate_simplex(d):
    return np.array([1/d] * d), np.array([[1] * d]), np.array([1]), d, 'simplex'

x, A, b, k, name = generate_simplex(5)
sparse_dikin = SparseDikinWalk(r = 0.9, thin = 1)
dikin_res = sparse_dikin.generateCompleteWalk(10_000, x, A, b, k, burn = 100)
```
We can also provide a corollary example using the dense formulation.
```python
import numpy as np
from scipy.sparse import csr_matrix, lil_matrix, csr_array
from polytopewalk.dense import DikinWalk
from polytopewalk import FacialReduction

def generate_simplex(d):
    return np.array([1/d] * d), np.array([[1] * d]), np.array([1]), d, 'simplex'

fr = FacialReduction()
_, A, b, k, name = generate_simplex(5)

polytope = fr.reduce(A, b, k, sparse = False)
dense_A = polytope.dense_A
dense_b = polytope.dense_b

dc = DenseCenter()
init = dc.getInitialPoint(dense_A, dense_b)

dikin_res = dikin.generateCompleteWalk(1_000, init, dense_A, dense_b, burn = 100)
```

## Testing
The `tests` folder includes comprehensives tests of the Facial Reduction algorithm, Initialization, Weights from MCMC algorithms, and Sparse/Dense Random Walk algorithms in both Python and C++. The user can run each of the files separately to make sure the package passes all of the test suites. 
