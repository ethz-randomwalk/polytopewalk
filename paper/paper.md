---
title: 'PolytopeWalk: Sparse MCMC Sampling over Polytopes'
tags:
  - MCMC methods
  - sparsity
  - interior-point methods
  - polytopes
  - facial reduction
authors:
  - name: Benny Sun
    affiliation: 1
  - name: Yuansi Chen
    affiliation: 2


affiliations:
 - name: Department of Statistics, Duke University
   index: 1
 - name: Department of Mathematics, ETH Zurich
   index: 2
date: 28 February 2025
bibliography: paper.bib
---

# Summary


High dimensional sampling is an important computational tool in statistics and other  disciplines, with applications ranging from Bayesian statistical uncertainty quantification, metabolic modeling in systems biology to volume computation. We present ``PolytopeWalk``, a new scalable Python library designed for uniform sampling over polytopes.  The library provides an end-to-end solution, that includes preprocessing algorithms such as facial reduction and initialization methods. Six state-of-the-art MCMC algorithms on polytopes are implemented, including the Dikin, Vaidya, and John Walk. Additionally, we introduce novel sparse constrained formulations of these algorithms, enabling efficient sampling from sparse polytopes of the form $\mathcal{K}_2 = \{x \in \mathbb{R}^d \ | \  Ax = b, x \succeq_k 0\}$. This implementation maintains sparsity in $A$, ensuring scalability to high dimensional settings $(d > 10^5)$. ``PolytopeWalk`` is  available at github.com/ethz-randomwalk/polytopewalk with documentation at polytopewalk.readthedocs.io/ .


# Acknowledgements

Much of the work was done while Yuansi Chen was an assistant professor in the Department of Statistical Science at Duke University. Both authors are partially supported by NSF CAREER Award DMS-2237322, Sloan Research Fellowship and Ralph E. Powe Junior Faculty Enhancement Awards. 

# References