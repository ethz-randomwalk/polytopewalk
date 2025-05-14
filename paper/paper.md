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
date: 4 March 2025
bibliography: paper.bib
---

# Summary

High dimensional sampling is an important computational tool in statistics, with applications in stochastic simulation, volume computation, and fast randomized algorithms. We present ``PolytopeWalk``, a scalable library designed for sampling from a uniform distribution over polytopes, which are bounded geometric objects formed by linear inequalities. For sampling, we use Markov chain Monte Carlo (MCMC) methods, defined as a family of algorithms for generating approximate samples from a target probability distribution. Six state-of-the-art MCMC algorithms are implemented, including the Dikin, Vaidya, and John Walk. Additionally, we introduce novel sparse constrained formulations of these algorithms, enabling efficient sampling from sparse polytopes of the form $\mathcal{K}_2 = \{x \in \mathbb{R}^d \ | \ Ax = b, x \succeq_k 0\}$. This implementation maintains sparsity in $A$, ensuring scalability to high dimensional settings $(d > 10^4)$. Finally, ``PolytopeWalk`` includes novel implementations of preprocessing algorithms such as facial reduction and initialization, thus providing an end-to-end solution.

# Statement of Need

High dimensional sampling is a fundamental problem in many computational disciplines such as statistics, probability, and operation research. For example, sampling is applied in portfolio optimization [@DBLP:journals/corr/abs-1803-05861], metabolic networks in biology [@COBRA] and volume approximation over convex shapes [@Simonovits2003]. Markov chain Monte Carlo (MCMC) sampling algorithms offer a natural and scalable solution to this problem. These algorithms construct a Markov chain whose stationary distribution matches the target distribution. By running the chain for a large number of steps to ensure mixing, MCMC algorithms can efficiently generate approximately independent samples close to the target distribution, while not suffering from the curse of dimension issues.

This package focuses on sampling from a uniform distribution over a user-specified polytope. We define the polytope as the following. Let $A \in \mathbb{R}^{n \times d}$, $b \in \mathbb{R}^n$ and let $x \succeq_k y$ mean that the last $k$-coordinates of $x$ are greater than or equal to the corresponding coordinates of $y$, i.e., $\{x_{d-k+1} - y_{d-k+1} \ge 0, ... , x_{d} - y_{d} \ge 0\}$. Depending on whether we allow equality constraints, the sampling problem can be formalized in two forms:

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

Large polytopes with sparse constraints are common in many applications. The largest human metabolic network RECON3D is modeled as a $13543$-dimensional sparse polytope [@10.1093/nar/gkv1049]. Moreover, linear programming datasets from `NetLib` are naturally in the constrained form, where $A$ matrix is sparse. These applications motivate the need for MCMC algorithms that leverage $\mathcal{K}_2$ form. We implement novel interior-point-method-based MCMC algorithms optimized for large and sparse constrained polytopes. By exploiting sparsity, our algorithms scale well in per-iteration cost and sampling efficiency as a function of increasing dimension, enabling effective sampling from polytopes with dimensions exceeding $10^4$.

Interior-point-method-based MCMC sampling algorithms on a polytope are modifications of the Ball Walk [@vempala2005], incorporating key concepts from interior-point methods in optimization. These algorithms operate in two primary steps. First, the algorithm generates a proposal distribution whose covariance matrix is state-dependent and equal to the inverse of the Hessian matrix of a specified barrier function, capturing the local geometry of the polytope. Second, the algorithm employs the Metropolis-Hastings accept-reject step to ensure that its stationary distribution is uniform on the polytope [@Metropolis1953; @10.1093/biomet/57.1.97]. Using a state-dependent proposal distribution that adapts to the polytope's local geometry, these MCMC algorithms achieve an improved mixing rate.

In ``PolytopeWalk``, we implement 4 interior-point-method-based MCMC sampling algorithms in both the sparse-constrained and full-dimensional formulation. ``PolytopeWalk`` makes meaningful strides in the open-source development of MCMC, speeding up calculations for sparse high-dimensional sampling. Finally, we provide an an open-source implementation of the Facial Reduction algorithm, described in detail in the Preprocessing Algorithms section. 

# Package Overview

``PolytopeWalk`` is an open-source library written in C++ with Python wrapper code, providing accelerated MCMC sampling algorithms in both $\mathcal{K}_1$ and $\mathcal{K}_2$ formulation. The source code is written in C++ with `Eigen` for linear algebra [@eigenweb], `glpk` for linear programming [@glpk], and `pybind` for Python binding [@pybind11]. In Python, ``PolytopeWalk`` relies on both NumPy [@harris2020array] and SciPy [@2020SciPy-NMeth].

![Code Structure of Package](images/Code_Design.pdf){ width=80% }

## Random Walk Algorithms

Mixing times refer to the required number of steps to converge to stationary distribution. In each, $d$ refers to the dimension of the polytope and $n$ refers to the number of boundaries ($\mathcal{K}_1$ dimensions). In the first 2 walks, $R^2/r^2$ means where the convex body contains a ball of radius $r$ and is mostly contained in a ball of radius $R$.

| Name              | Mixing Time             | Author              |
|:----------------:|:---------------------:|:-------------------:|
| `Ball Walk`       | $O(d^2R^2/r^2)$     | [Vempala (2005)](https://faculty.cc.gatech.edu/~vempala/papers/survey.pdf) |
| `Hit and Run`     | $O(d^2R^2/r^2)$     | [Lovasz (1999)](https://link.springer.com/content/pdf/10.1007/s101070050099.pdf) |
| `Dikin Walk`      | $O(nd)$             | [Sachdeva et al. (2015)](https://arxiv.org/pdf/1508.01977) |
| `Vaidya Walk`     | $O(n^{1/2}d^{3/2})$ | [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html) |
| `John Walk`       | $O(d^{2.5})$        | [Chen et al. (2018)](https://jmlr.org/papers/v19/18-158.html) |
| `Lee Sidford Walk`| $O(d^{2})$         | [Laddha et al. (2019)](https://arxiv.org/abs/1911.05656) |

## Preprocessing Algorithms

``PolytopeWalk`` comes with 2 preprocessing algorithms: initialization and facial reduction.

**Initialization:** If the user cannot specify a point inside of the polytope to start, ``PolytopeWalk`` provides a class to compute an initial point well within the polytope for both the full-dimensional formulation and constrained formulation.

**Facial Reduction:** We adopt the facial reduction algorithm implementation from Drusvyatskiy's research [@drusvyatskiy2017many; @im2023revisiting]. In the constrained formulation $\mathcal{K}_2 = \{x \in \mathbb{R}^d \ | \ Ax = b, x \succeq_k 0\}$, degeneracy occurs when there is a lack of strict feasibility in the polytope: there does not exist an $x \in \mathbb{R}^d$ such that $Ax = b$ and $x \succ_k 0$. Thus, degeneracy exists in polytopes when the lower-dimensional polytope is embedded in a higher dimension. The facial reduction algorithm eliminates variables in the last k dimensions fixed at $0$, thus ensuring numerical stability for sampling. 

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


Table II contrasts the features of ``PolytopeWalk`` with `Volesti` [@Chalkis_2021], `WalkR` [@Yao2017], and `Polyrun` [@CIOMEK2021100659]. `Volesti` is implemented in C++ with some of its code represented in the Python library `Dingo`. `Polyrun` only works on Java and `WalkR` on R. Notably, `WalkR` was removed from the CRAN repository, motivating further open source MCMC sampling development. 

# Acknowledgements

Much work was done while Yuansi Chen was an assistant professor in the Department of Statistical Science at Duke University. Both authors are partially supported by NSF CAREER Award DMS-2237322, Sloan Research Fellowship and Ralph E. Powe Junior Faculty Enhancement Awards. 

# References