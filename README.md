![example branch parameter](https://github.com/ethz-randomwalk/polytopewalk/actions/workflows/ciwheels.yml/badge.svg?branch=main)
# PolytopeWalk
**PolytopeWalk** is a `C++` library for running MCMC sampling algorithms to generate samples from a uniform distribution over a polytope with a `Python` interface. It handles preprocessing of the polytope and initialization as well. Current implementations include the Dikin Walk, John Walk, Vaidya Walk, Ball Walk, Lee Sidford Walk, and Hit-and-Run in both the full-dimensional formulation and the sparse constrained formulation. Code includes facial reduction and initialization algorithms for pre-processing as well. Sample code that samples from both real polytopes from a data set and artificial polytopes are shown in the Examples folder.

## Developer Installation Instructions 
First, we need to install package prerequisites (listed in each of the operating systems)
- macOS: ``brew install eigen glpk``
- Windows: ``choco install glpk eigen``
- Linux: ``yum install -y epel-release eigen3-devel glpk-devel``

Finally, we can install **PolytopeWalk** via pip:
```
git clone https://github.com/ethz-randomwalk/polytopewalk
cd polytopewalk
pip install .
```
