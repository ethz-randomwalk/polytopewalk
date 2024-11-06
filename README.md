![example branch parameter](https://github.com/ethz-randomwalk/polytopewalk/actions/workflows/ciwheels.yml/badge.svg?branch=main)
# PolytopeWalk
**PolytopeWalk** is a `C++` library for running MCMC sampling algorithms to generate samples from a uniform distribution over a polytope with a `Python` interface. It handles preprocessing of the polytope and initialization as well. Current implementations include the Dikin Walk, John Walk, Vaidya Walk, Ball Walk, Lee Sidford Walk, and Hit-and-Run in both the full-dimensional formulation and the sparse constrained formulation. Code includes facial reduction and initialization algorithms for pre-processing as well. Sample code that samples from both real polytopes from a data set and artificial polytopes are shown in the Examples folder.

# Installation

## Dependencies
polytopewalk requires:
- Python (>= 3.9)
- NumPy (>= 1.20)
- SciPy (>= 1.6.0)

## User installation
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


### Complie C++ from source (not necessary)
Only do this, if there is need to run and test C++ code directly. For normal users, we recommend only using the Python interface. 

Build with cmake
```bash
git clone https://github.com/ethz-randomwalk/polytopewalk.git && cd polytopewalk
cmake -B build -S .
make
sudo make install
```