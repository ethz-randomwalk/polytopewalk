![example branch parameter](https://github.com/ethz-randomwalk/polytopewalk/actions/workflows/ciwheels.yml/badge.svg?branch=main)
# PolytopeWalk
**PolytopeWalk** is a `C++` library for running MCMC sampling algorithms to generate samples from a uniform distribution over a polytope with a `Python` interface. It handles preprocessing of the polytope and initialization as well. Current implementations include the Dikin Walk, John Walk, Vaidya Walk, Ball Walk, Lee Sidford Walk, and Hit-and-Run. Sample code that samples from both real polytopes from a data set and artificial polytopes are shown in the Example folder.

## Developer Installation Instructions 
First, we need to install package prerequisites (listed in each of the operating systems)
- macOS: ``brew install eigen ipopt``
- Windows: ``vcpkg install eigen3 coin-or-ipopt``
- Linux: ``yum install -y epel-release eigen3-devel coin-or-Ipopt-devel``

Next, we need to install ifopt from its source on Github: 
```
git clone https://github.com/ethz-adrl/ifopt.git && cd ifopt
mkdir build && cd build
cmake ..
make
sudo make install
```

Finally, we can install **PolytopeWalk** via pip:
```
git clone https://github.com/bsun1220/polytopewalk
cd polytopewalk
pip install .
```
