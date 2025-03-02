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

High dimensional sampling is an important computational tool in statistics and other  disciplines, with applications ranging from Bayesian statistical uncertainty quantification, metabolic modeling in systems biology, to volume computation. We present ``PolytopeWalk``, a new scalable Python library designed for uniform sampling over polytopes.  The library provides an end-to-end solution, that includes preprocessing algorithms such as facial reduction and initialization methods. Six state-of-the-art MCMC algorithms on polytopes are implemented, including the Dikin, Vaidya, and John Walk. Additionally, we introduce novel sparse constrained formulations of these algorithms, enabling efficient sampling from sparse polytopes of the form $\mathcal{K}_2 = \{x \in \mathbb{R}^d \ | \  Ax = b, x \succeq_k 0\}$. This implementation maintains sparsity in $A$, ensuring scalability to high dimensional settings $(d > 10^4)$. We also include novel implementation of several pre-processing steps, including the Facial Reduction algorithm. ``PolytopeWalk`` is  available at [github.com/ethz-randomwalk/polytopewalk](github.com/ethz-randomwalk/polytopewalk) with documentation at [polytopewalk.readthedocs.io/](polytopewalk.readthedocs.io/).

# Statement of Need

High dimensional sampling is a fundamental problem in many computational disciplines such as statistics, probability, and operation research. For example, sampling is applied in portfolio optimization (@DBLP:journals/corr/abs-1803-05861), metabolic networks in systems biology (@COBRA) and volume approximation over convex shapes (@Simonovits2003). Markov chain Monte Carlo (MCMC) sampling algorithms offer a natural and scalable solution to this problem. These algorithms construct a Markov chain whose stationary distribution matches the target distribution. By running the chain for a large number of steps to ensure mixing, MCMC algorithms can efficiently generate approximately independent samples close to the target distribution, while not suffering from the curse of dimension issues.

This package focuses on sampling from a uniform distribution over a user-specified polytope. We define the polytope as the following. Let $A \in \mathbb{R}^{n \times d}$, $b \in \mathbb{R}^n$ and define  $x \succeq_k y$ to mean that the last $k$-coordinates of $x$ are greater than or equal to the corresponding coordinates of $y$, i.e., $\{x_{d-k+1} - y_{d-k+1} \ge 0, ... , x_{d} - y_{d} \ge 0\}$. Depending on whether we allow equality constraints, the sampling problem can be formalized in two forms:

\begin{enumerate}
    \item The full-dimensional form:
    \begin{align}
        \mathcal{K}_1 = \{x \in \mathbb{R}^d \ | Ax \le b\},
        \label{eq:full_dim}
    \end{align}
    where $\mathcal{K}_1$ is specified via $n$ inequality constraints. 
    \item The constrained form:
    \begin{align}
        \mathcal{K}_2 = \{x \in \mathbb{R}^d \ | \ Ax = b, x \succeq_k 0\},
        \label{eq:constrained}
    \end{align}
    where $\mathcal{K}_2$ is specified via $n$ equality constraints and $k$ coordinate inequality constraints. 
\end{enumerate}

Large polytopes with sparse constraints are common in many applications. For instance, the largest human metabolic network RECON3D is modeled as a $13543$-dimensional sparse polytope (@10.1093/nar/gkv1049). Additionally, linear programming datasets from `NetLib` are naturally in the constrained form, where $A$ matrix is sparse. These applications motivate the need for MCMC algorithms that leverage the large and sparse $\mathcal{K}_2$ formulation. We develop novel interior-point-method-based MCMC algorithms optimized for large and sparse constrained polytopes. By exploiting sparsity, our algorithms scale well in both per-iteration cost and sampling efficiency as a function of increasing dimension, enabling effective sampling from polytopes with dimensions exceeding $10^4$.

Interior-point-method-based MCMC sampling algorithms on a polytope are modifications of the Ball Walk (@vempala2005), incorporating key concepts from interior-point methods in optimization. These algorithms operate in two primary steps. First, the algorithm generates a proposal distribution whose covariance matrix is state-dependent and equal to the inverse of the Hessian matrix of a specified barrier function, capturing the local geometry of the polytope. Second, the algorithm employs the Metropolis-Hastings accept-reject step to ensure that its stationary distribution is uniform on the polytope. Using a state-dependent proposal distribution that adapts to the polytope's local geometry, these MCMC algorithms achieve an improved mixing rate.  Specific algorithms in this class include the Dikin Walk (@DBLP:journals/corr/SachdevaV15), Vaidya Walk (@JMLR:v19:18-158), John Walk (@JMLR:v19:18-158), and Lee Sidford Walk (@DBLP:journals/corr/abs-1911-05656), with theoretical guarantees therein.  Each of these methods leverages different barrier functions and are specialized for sampling distributioned truncated on a polytope.

In ``PolytopeWalk``, we implement 4 interior-point-method-based MCMC sampling algorithms in both the sparse, constrained formulation and the full-dimensional form. ``PolytopeWalk`` makes meaningful strides in the open-source development of MCMC, speeding up calculations for sparse high-dimensional sampling. On our Github, we demonstrate the scalability of ``PolytopeWalk`` against packages like `Volesti`. Finally, we provide an an open-source implementation of the Facial Reduction algorithm, used to handle degeneracy in polytopes.

# Package Overview

``PolytopeWalk`` is an open-source library written in C++ with Python wrapper code. There are 3 main components of the package: preprocessing, sampling, and post processing. ```PolytopeWalk``` provides accelerated MCMC sampling algorithms in both the $\mathcal{K}_1$ formulation~\eqref{eq:full_dim} and $\mathcal{K}_2$ formulation~\eqref{eq:constrained}. The source code is written in low-level C++ with `Eigen` for optimized linear algebra operations, `glpk` for a fast sparse linear programming solver, and `pybind` to enable Python binding. 

![Code Structure of Package](images/Code_Design.pdf){ width=100% }

## Random Walk Algorithms

There are 2 main submodules: `dense` and `sparse`. The `dense` module includes implementations of the 6 random walks algorithms in the full-dimensional formulation. The `sparse` module includes implementations in the sparse, constrained formulation. Table I provides a broad overview of their respective mixing rates and cited paper. 

| Name              | Mixing Time             | Author              |
|:----------------:|:---------------------:|:-------------------:|
| `Ball Walk`       | $\tau(d^2R^2/r^2)$     | [Vempala (2005)](https://faculty.cc.gatech.edu/~vempala/papers/survey.pdf) |
| `Hit and Run`     | $\tau(d^2R^2/r^2)$     | [Lovasz (1999)](https://link.springer.com/content/pdf/10.1007/s101070050099.pdf) |
| `Dikin Walk`      | $\tau(nd)$             | [Sachdeva et al. (2015)](https://arxiv.org/pdf/1508.01977) |
| `Vaidya Walk`     | $\tau(n^{1/2}d^{3/2})$ | [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html) |
| `John Walk`       | $\tau(d^{2.5})$        | [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html) |
| `Lee Sidford Walk`| $\tau(d^{2})$         | [Laddha et al. (2019)](https://arxiv.org/abs/1911.05656) |

The mixing times refer to the required number of steps to converge to stationary distribution. In each, $d$ refers to the dimension of the polytope and $n$ refers to the number of boundaries ($\mathcal{K}_1$ dimensions). In the first 2 walks, $R^2/r^2$ means where the convex body contains a ball of radius $r$ and is mostly contained in a ball of radius $R$.

## Preprocessing Algorithms

``PolytopeWalk`` comes with 2 preprocessing algorithms: initialization and facial reduction.

**Initialization:** If the user cannot specify a point inside of the polytope to start, ``PolytopeWalk`` provides a class to compute an initial point well within the polytope for both the full-dimensional formulation and constrained formulation. In $\mathcal{K}_2$, this is akin to solving the linear program maximizing $\delta \in \mathbb{R}$ s.t. $Ax =b$ and $x \succeq_k \vec{\delta}$. In $\mathcal{K}_1$, this is equivalent to maximizing $\delta \in \mathbb{R}$ such that $Ax + \delta \cdot \mathbb{1} \le b$ where $\mathbb{1}$ is the all-ones vector $\in \mathbb{R}^n$.

**Facial Reduction:** We adopt the facial reduction algorithm implementation from Drusvyatskiy's research (@drusvyatskiy2017many) (@im2023revisiting). The facial reduction algorithm deals with cases of degeneracy in the polytope. In the constrained formulation $\mathcal{K}_2 = \{x \in \mathbb{R}^d \ | \ Ax = b, x \succeq_k 0\}$, degeneracy occurs when there is a lack of strict feasibility in the polytope: there does not exist an $x \in \mathbb{R}^d$ such that $Ax = b$ and $x \succ_k 0$. The facial reduction algorithm eliminates variables in the last k dimensions fixed at $0$. It is far more commonly used in linear programming than MCMC sampling. Degeneracy exists in polytopes when the lower-dimensional polytope is embedded in a higher dimension. For example, if one attempts to use MCMC algorithms to sample from a 2-dimensional hypercube in 5 dimensions, the proposal distributions which generate an ellipsoid in 5 dimensions will certainly reject each proposal. It is crucial that one projects the polytope back into the 2 dimensional form to run the MCMC sampler efficiently. Facial Reduction accomplishes this task, ensuring numerical stability for sampling.

## Package Comparison

| Feature                  | ``PolytopeWalk`` | `Volesti` | `WalkR` | `Polyrun` |
|--------------------------|:-------------:|:--------:|:------:|:--------:|
| Constrained Formulation  | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{green}Y}$$ | $${\color{green}Y}$$ |
| Sparse Friendly          | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ |
| C++ Implementation       | $${\color{green}Y}$$ | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ |
| Facial Reduction         | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ |
| Dikin Walk               | $${\color{green}Y}$$ | $${\color{green}Y}$$ | $${\color{green}Y}$$ | $${\color{red}N}$$ |
| Vaidya Walk              | $${\color{green}Y}$$ | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ |
| John Walk                | $${\color{green}Y}$$ | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ |
| Lee-Sidford Walk         | $${\color{green}Y}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ | $${\color{red}N}$$ |


Table II contrasts the features of ``PolytopeWalk`` with `Volesti` (@Chalkis_2021), `WalkR` (@Yao2017), and `Polyrun` (@CIOMEK2021100659). ``PolytopeWalk`` is the first open-source package that enables users to leverage sparsity and use constrained formulations of the polytopes to generate points uniformly with fast MCMC algorithms. We are also one of the first to primarily focus on barrier walk MCMC samplers. ``PolytopeWalk`` includes a C++ implementation with corresponding Python wrapper code. Conversely, `Volesti` is implemented in C++ with some of its code represented in the Python library `Dingo`. `Polyrun` only works on Java and `WalkR` on R. Notably, `WalkR` was removed from the CRAN repository, motivating further open source MCMC sampling development. 

# Acknowledgements

Much of the work was done while Yuansi Chen was an assistant professor in the Department of Statistical Science at Duke University. Both authors are partially supported by NSF CAREER Award DMS-2237322, Sloan Research Fellowship and Ralph E. Powe Junior Faculty Enhancement Awards. 

# References