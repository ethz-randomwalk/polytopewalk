#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "utils/FullWalkRun.hpp"
namespace py = pybind11;

template <class RandomWalkBase = RandomWalk> class PyRandomWalk : public RandomWalkBase {
public:
    using RandomWalkBase::RandomWalkBase; // Inherit constructors
    MatrixXd generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burnin, int thin, int seed) override{
            PYBIND11_OVERRIDE_PURE(
                MatrixXd,
                RandomWalkBase,
                generateCompleteWalk,
                niter,
                init,
                A,
                b,
                burnin,
                thin, 
                seed
            );
    }
};

template <class BarrierWalkBase = BarrierWalk> class PyBarrierWalk: public PyRandomWalk<BarrierWalkBase> {
public:
    using PyRandomWalk<BarrierWalkBase>::PyRandomWalk;
    MatrixXd generateCompleteWalk(const int niter, VectorXd& init, const MatrixXd& A, const VectorXd& b, int burnin, int thin, int seed) override{
            PYBIND11_OVERRIDE(
                MatrixXd,
                BarrierWalkBase,
                generateCompleteWalk,
                niter,
                init,
                A,
                b,
                burnin,
                thin, 
                seed
            );
    }
    VectorXd generateWeight(const VectorXd& x, const MatrixXd& A, const VectorXd& b) override{
            PYBIND11_OVERRIDE(
                VectorXd,
                BarrierWalkBase,
                generateWeight,
                x,
                A,
                b
            );
        }
    void setDistTerm(int d, int n) override{
            PYBIND11_OVERRIDE(
                void,
                BarrierWalkBase,
                setDistTerm,
                d,
                n
            );
    }
};


template <class SparseRandomWalkBase = SparseRandomWalk> class PySparseRandomWalk : public SparseRandomWalkBase {
public:
    using SparseRandomWalkBase::SparseRandomWalkBase; // Inherit constructors
    MatrixXd generateCompleteWalk(
        const int niter, 
        const VectorXd& init, 
        const SparseMatrixXd& A, 
        const VectorXd& b, 
        int k, 
        int burnin,
        int thin, 
        int seed
        ) override
    {
            PYBIND11_OVERRIDE_PURE(
                MatrixXd,
                SparseRandomWalkBase,
                generateCompleteWalk,
                niter,
                init,
                A,
                b,
                k,
                burnin,
                thin, 
                seed
            );
    }
};
template <class SparseBarrierWalkBase = SparseBarrierWalk> class PySparseBarrierWalk: public PySparseRandomWalk<SparseBarrierWalkBase> {
public:
    using PySparseRandomWalk<SparseBarrierWalkBase>::PySparseRandomWalk; 
    MatrixXd generateCompleteWalk(
        const int niter, 
        const VectorXd& init, 
        const SparseMatrixXd& A, 
        const VectorXd& b,
        int k,
        int burnin,
        int thin, 
        int seed
        ) override
        {
            PYBIND11_OVERRIDE(
                MatrixXd,
                SparseBarrierWalkBase,
                generateCompleteWalk,
                niter,
                init,
                A,
                b,
                k,
                burnin,
                thin, 
                seed
            );
    }

    VectorXd generateWeight(const VectorXd& x, const SparseMatrixXd& A, int k) override{
            PYBIND11_OVERRIDE(
                VectorXd,
                SparseBarrierWalkBase,
                generateWeight,
                x,
                A,
                k
            );
    }
    void setDistTerm(int d, int n) override{
            PYBIND11_OVERRIDE(
                void,
                SparseBarrierWalkBase,
                setDistTerm,
                d,
                n
            );
    }
};


PYBIND11_MODULE(polytopewalk, m) {
    m.doc() = "Polytopewalk Library";

    m.def("denseFullWalkRun", &denseFullWalkRun, 
    R"doc(
    Dense Central Function. Starts with polytope in sparse, constrained formulation. Computes facial reduction to preprocess polytope and uses the dense, full-dimensional polytope to run MCMC sampler before converting it back into original formulation. 

    Parameters
    ----------
    niter : int
        Number of iterations.
    A : numpy.ndarray
        Constraint matrix.
    b : numpy.ndarray
        Constraint vector.
    k : int
        Dimensionality of the polytope.
    walk : RandomWalk
        Random walk instance.
    fr : FacialReduction
        Facial reduction object.
    dc : DenseCenter
        Dense center object.
    burnin : int, optional
        Number of burn-in steps (default is 0).
    thin : int, optional
        Number of samples to thin (default is 1).
    seed : int, optional
        Seed number for reproducibility (default is -1 meaning no fixed setting).

    Returns
    --------
    numpy.ndarray
        List of sampled points.
    )doc",
    py::arg("niter"), py::arg("A"), py::arg("b"), py::arg("k"), py::arg("walk"), py::arg("fr"), py::arg("dc"), py::arg("burnin") = 0, py::arg("thin") = 1, py::arg("seed") = -1);

    m.def("sparseFullWalkRun", &sparseFullWalkRun, 
    R"doc(
    Sparse Central Function. Starts with polytope in sparse, constrained formulation. Computes facial reduction to preprocess polytope and uses the reduced constrained polytope to run MCMC sampler before converting it back into original formulation.

    Parameters
    ----------
    niter : int
        Number of iterations.
    A : numpy.ndarray
        Constraint matrix.
    b : numpy.ndarray
        Constraint vector.
    k : int
        Dimensionality of the polytope.
    walk : RandomWalk
        Random walk instance.
    fr : FacialReduction
        Facial reduction object.
    sc : SparseCenter
        Sparse center object.
    burnin : int, optional
        Number of burn-in steps (default is 0).
    thin : int, optional
        Number of samples to thin (default is 1).
    seed : int, optional
        Seed number for reproducibility (default is -1 meaning no fixed setting).

    Returns
    --------
    numpy.ndarray
        List of sampled points.
    )doc", 
    py::arg("niter"), py::arg("A"), py::arg("b"), py::arg("k"), py::arg("walk"), py::arg("fr"), py::arg("sc"), py::arg("burnin") = 0, py::arg("thin") = 1, py::arg("seed") = -1);

    auto m_dense = m.def_submodule("dense", "Dense Module");
    auto m_sparse = m.def_submodule("sparse", "Sparse Module");

    py::class_<DenseCenter>(m_dense, "DenseCenter", "Initialization Algorithm for Dense Polytopes.")
        .def(py::init<>(), "Initialization for Center Algorithm.")
        .def("getInitialPoint", &DenseCenter::getInitialPoint, 
            R"doc(
            Finds analytical center for Ax <= b.

            Parameters
            ----------
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.

            Returns
            --------
            numpy.ndarray
                Point well within polytope.
            )doc",
            py::arg("A"), py::arg("b"));
    
    py::class_<SparseCenter>(m_sparse, "SparseCenter", "Initialization Algorithm for Sparse Polytopes.")
        .def(py::init<>(), "Initialization for Center Algorithm.")
        .def("getInitialPoint",&SparseCenter::getInitialPoint,
            R"doc(
            Finds analytical center Ax = b, x >=_k 0.

            Parameters
            -----------
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.
            k : int
                Dimensionality of the polytope.

            Returns
            --------
            numpy.ndarray
                Point well within polytope.
            )doc",
            py::arg("A"), py::arg("b"), py::arg("k"));
    
    py::class_<RandomWalk, PyRandomWalk<>>(m_dense, "RandomWalk", "Random Walk Superclass Implementation.")
        .def(py::init<>(), 
            R"doc(
            Initialization for Random Walk Super Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.
            
            )doc")
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk, 
            R"doc(
            Generate values from Random Walk (virtual function).

            Parameters
            -----------
            niter : int
                Number of iterations.
            init : numpy.ndarray
                Initial point to start sampling from.
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.
            burnin : int, optional
                Constant for how many to exclude initially (default is 0).
            thin : int, optional
                Number of samples to thin (default is 1).
            seed : int, optional
                Seed number for reproducibility (default is -1 meaning no fixed setting).

            Returns
            --------
            numpy.ndarray
                List of sampled points.
            )doc", 
            py::arg("niter"), py::arg("init"), py::arg("A"), py::arg("b"), py::arg("burnin") = 0, py::arg("thin") = 1, py::arg("seed") = -1
        );
    
    py::class_<BallWalk, RandomWalk>(m_dense, "BallWalk", "Ball Walk Implementation.")
        .def(py::init<double>(), 
            R"doc(
            Initialization for Ball Walk Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for ball (default is 0.5).

            )doc",
            py::arg("r") = 0.5);
    
    py::class_<HitAndRun, RandomWalk>(m_dense, "HitAndRun", "Hit-Run Implementation.")
        .def(py::init<double, double>(),  
            R"doc(
            Initialization for Hit and Run Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for starting distance (default is 0.5).
            err : double, optional
                Constant for closeness to edge of polytope (default is 0.01).

            )doc",
            py::arg("r") = 0.5, py::arg("err") = 0.01);

    py::class_<BarrierWalk, RandomWalk, PyBarrierWalk<>>(m_dense, "BarrierWalk", "Barrier Walk Implementation.")
        .def(py::init<double>(), 
            R"doc(
            Initialization for Barrier Walk Super Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for starting distance (default is 0.5).

            )doc",
            py::arg("r") = 0.5)
        .def("generateWeight", &BarrierWalk::generateWeight, 
            R"doc(
            Generate weight from Barrier Walk (virtual function).

            Parameters
            -----------
            x : numpy.ndarray
                Point inside polytope.
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.

            Returns
            --------
            numpy.ndarray
                Weight vector (specified by walk type). 
            )doc", 
            py::arg("x"), py::arg("A"), py::arg("b"))
        .def("generateCompleteWalk", &RandomWalk::generateCompleteWalk, 
            R"doc(
            Generate values from Barrier Walk (virtual function).

            Parameters
            -----------
            niter : int
                Number of steps to sample from.
            init : numpy.ndarray
                Initial point to start sampling from.
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.
            burnin : int, optional
                Constant for how many to exclude initially (default is 0).
            thin : int, optional
                Number of samples to thin (default is 1).
            seed : int, optional
                Seed number for reproducibility (default is -1 meaning no fixed setting).

            Returns
            --------
            numpy.ndarray
                List of sampled points.
            )doc", 
            py::arg("niter"), py::arg("init"), py::arg("A"), py::arg("b"), py::arg("burnin") = 0, py::arg("thin") = 1, py::arg("seed") = -1
        );
    
    py::class_<DikinWalk, BarrierWalk, PyBarrierWalk<DikinWalk>>(m_dense, "DikinWalk", "Dikin Walk Implementation.")
        .def(py::init<double>(), 
            R"doc(
            Initialization for Dikin Walk Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for Dikin Ellipsoid (default is 0.5).
            )doc",
            py::arg("r") = 0.5);
    
    py::class_<VaidyaWalk, BarrierWalk, PyBarrierWalk<VaidyaWalk>>(m_dense, "VaidyaWalk", "Vaidya Walk Implementation.")
       .def(py::init<double>(), 
            R"doc(
            Initialization for Vaidya Walk Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for Vaidya Ellipsoid (default is 0.5).
            )doc",
            py::arg("r") = 0.5);
    
    py::class_<DikinLSWalk, BarrierWalk, PyBarrierWalk<DikinLSWalk>>(m_dense, "DikinLSWalk", "Lee Sidford Walk Implementation.")
        .def(py::init<double, double, double, int>(), 
            R"doc(
            Initialization for Lee Sidford Walk Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for Lee-Sidford Ellipsoid (default is 0.5).
            g_lim : double, optional
                Constant for stopping gradient norm in gradient descent (default is 0.01).
            step_size : double, optional
                Constant for step size in gradient descent (default is 0.1).
            max_iter : int, optional
                Constant for maximum number of gradient descent iterations (default is 1000).
            )doc",
            py::arg("r") = 0.5, py::arg("g_lim") = 0.01, py::arg("step_size") = 0.1, 
            py::arg("max_iter") = 1000);
    
    py::class_<JohnWalk, BarrierWalk, PyBarrierWalk<JohnWalk>>(m_dense, "JohnWalk", "John Walk Implementation.")
        .def(py::init<double, double, int>(), 
            R"doc(
            Initialization for John Walk Class.
            Runs on Full-Dimensional Polytope Form: Ax <= b.

            Parameters
            -----------
            r : double, optional
                Radius for John Ellipsoid (default is 0.5).
            lim : double, optional
                Constant for stopping limit in fixed-point iteration (default is 1e-5).
            max_iter : int, optional
                Constant for maximum number of fixed point iterations (default is 1000).
            )doc",
            py::arg("r") = 0.5, py::arg("lim") = 1e-5, py::arg("max_iter") = 1000);
    
    py::class_<FacialReduction>(m, "FacialReduction", "Facial Reduction Implementation.")
        .def(py::init<double>(), 
            R"doc(
            Initialization for Facial Reduction Class.

            Parameters
            -----------
            err_dc : double, optional
                Error sensitivity for decomposition calculation (default is 1e-6).
            )doc",
            py::arg("err_dc") = 1e-6)
        .def("reduce", &FacialReduction::reduce, 
            R"doc(
            Completes facial reduction on Ax = b, x >=_k 0. 

            Parameters
            -----------
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.
            k : int
                Dimensionality of polytope.
            sparse : bool
                Includes only sparse constrained polytope or adds dense full-dimensional polytope.

            Returns
            --------
            FROutput
                Facial Reduction results object.
            )doc",
            py::arg("A"), py::arg("b"), py::arg("k"), py::arg("sparse"));
    
    py::class_<FROutput>(m, "FROutput", "Output for Facial Reduction.")
        .def_readwrite("sparse_A", &FROutput::sparse_A, "Constrained form Ax = b, x >=_k 0.")
        .def_readwrite("sparse_b", &FROutput::sparse_b, "Constrained form Ax = b, x >=_k 0.")
        .def_readwrite("saved_V", &FROutput::saved_V, "PAVv = Pb decomposition.")
        .def_readwrite("dense_A", &FROutput::dense_A, "Full-dim form Ax <= b.")
        .def_readwrite("dense_b", &FROutput::dense_b, "Full-dim form Ax <= b.")
        .def_readwrite("Q", &FROutput::Q, "Matrix used to go between forms.")
        .def_readwrite("z1", &FROutput::z1, "Vector used to go between forms.");
    
    py::class_<SparseRandomWalk, PySparseRandomWalk<>>(m_sparse, "SparseRandomWalk", "Sparse Random Walk Super Class Implementation.")
        .def(py::init<double>(), 
            R"doc(
            Initialization for Sparse Random Walk Super Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            err : double, optional
                Constant for error term term (default is 1e-6).
            )doc",
            py::arg("err") = 1e-6)
            .def("generateCompleteWalk", &SparseRandomWalk::generateCompleteWalk, 
            R"doc(
            Generate values from Sparse Random Walk (virtual function).

            Parameters
            -----------
            niter : int
                Number of steps to sample from.
            init : numpy.ndarray
                Initial point to start sampling from.
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.
            k : int
                Dimensionality of polytope.
            burnin : int, optional
                Constant for how many to exclude initially (default is 0).
            thin : int, optional
                Number of samples to thin (default is 1).
            seed : int, optional
                Seed number for reproducibility (default is -1 meaning no fixed setting).

            Returns
            --------
            numpy.ndarray
                List of sampled points.
            )doc",
            py::arg("niter"), py::arg("init"), py::arg("A"), py::arg("b"), py::arg("k"), py::arg("burnin") = 0, py::arg("thin") = 1, py::arg("seed") = -1
            );
    
    py::class_<SparseBallWalk, SparseRandomWalk>(m_sparse, "SparseBallWalk", "Sparse Ball Walk Implementation.")
        .def(py::init<double>(), 
            R"doc(
            Initialization for Sparse Ball Walk Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            r : double, optional
                Radius for ball (default is 0.5).

            )doc",
            py::arg("r") = 0.5);
    
    py::class_<SparseHitAndRun, SparseRandomWalk>(m_sparse, "SparseHitAndRun", "Sparse Hit and Run Implementation.")
        .def(py::init<double, double>(),  
            R"doc(
            Initialization for Sparse Hit and Run Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            r : double, optional
                Radius for starting distance (default is 0.5).
            err : double, optional
                Constant for closeness to edge of polytope (default is 0.01).

            )doc",
            py::arg("r") = 0.5, py::arg("err") = 0.01);

    py::class_<SparseBarrierWalk, SparseRandomWalk, PySparseBarrierWalk<>>(m_sparse, "SparseBarrierWalk", "Sparse Barrier Walk Implementation.")
        .def(py::init<double, double>(), 
            R"doc(
            Initialization for Sparse Barrier Walk Super Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            r : double, optional
                Radius for starting distance (default is 0.5).
            err : double, optional
                Constant for error term in g^{-1}(x) (default is 1e-6).

            )doc",
            py::arg("r") = 0.5, py::arg("err") = 1e-6)
        .def("generateWeight", &SparseBarrierWalk::generateWeight, 
            R"doc(
            Generate weight from Sparse Barrier Walk (virtual function).

            Parameters
            ----------
            x : numpy.ndarray
                Point inside polytope.
            A : numpy.ndarray
                Constraint matrix.
            k : int
                Dimensionality of polytope.

            Returns
            --------
            numpy.ndarray
                Weight vector (specified by walk type). 
            )doc",
            py::arg("x"), py::arg("A"), py::arg("k"))
        .def("generateCompleteWalk", &SparseRandomWalk::generateCompleteWalk, 
            R"doc(
            Generate values from Sparse Barrier Walk.

            Parameters
            -----------
            niter : int
                Number of steps to sample from.
            init : numpy.ndarray
                Initial point to start sampling from.
            A : numpy.ndarray
                Constraint matrix.
            b : numpy.ndarray
                Constraint vector.
            k : int
                Dimensionality of polytope.
            burnin : int, optional
                Constant for how many to exclude initially (default is 0).
            thin : int, optional
                Number of samples to thin (default is 1).
            seed : int, optional
                Seed number for reproducibility (default is -1 meaning no fixed setting).

            Returns
            --------
            numpy.ndarray
                List of sampled points.
            )doc",  
            py::arg("niter"), py::arg("init"), py::arg("A"), py::arg("b"), py::arg("k"), py::arg("burnin") = 0, py::arg("thin") = 1, py::arg("seed") = -1
        );
    
    py::class_<SparseDikinWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseDikinWalk>>(m_sparse, "SparseDikinWalk", "Sparse Dikin Walk Implementation.")
        .def(py::init<double, double>(), 
            R"doc(
            Initialization for Sparse Dikin Walk Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            r : double, optional
                Radius for Dikin Ellipsoid (default is 0.5).
            err : double, optional
                Constant for error term in g^{-1}(x) (default is 1e-6).
            )doc",
            py::arg("r") = 0.5, py::arg("err") = 1e-6);
    
    py::class_<SparseVaidyaWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseVaidyaWalk>>(m_sparse, "SparseVaidyaWalk", "Sparse Vaidya Walk Implementation.")
        .def(py::init<double, double>(), 
            R"doc(
            Initialization for Sparse Vaidya Walk Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            r : double, optional
                Radius for Vaidya Ellipsoid (default is 0.5).
            err : double, optional
                Constant for error term in g^{-1}(x) (default is 1e-6).
            )doc",
             py::arg("r") = 0.5, py::arg("err") = 1e-6);
    
    py::class_<SparseJohnWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseJohnWalk>>(m_sparse, "SparseJohnWalk", "Sparse John Walk Implementation.")
        .def(py::init<double, double, int, double>(), 
            R"doc(
            Initialization for Sparse John Walk Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.

            Parameters
            -----------
            r : double, optional
                Radius for John Ellipsoid (default is 0.5).
            lim : double, optional
                Constant for stopping limit in fixed-point iteration (default is 1e-5).
            max_iter : int, optional
                Constant for maximum number of fixed point iterations (default is 1000).
            err : double, optional
                Constant for error term in g^{-1}(x) (default is 1e-6).
            )doc", 
            py::arg("r") = 0.5, py::arg("lim") = 1e-5, 
            py::arg("max_iter") = 1000, py::arg("err") = 1e-6);
    
    py::class_<SparseDikinLSWalk, SparseBarrierWalk, PySparseBarrierWalk<SparseDikinLSWalk>>(m_sparse, "SparseDikinLSWalk", "Sparse Lee Sidford Walk Implementation.")
        .def(py::init<double, double, double, int, double>(), 
            R"doc(
            Initialization for Sparse Lee Sidford Walk Class.
            Runs on Constrained Polytope Form: Ax = b, x >=_k 0.
            
            Parameters
            -----------
            r : double, optional
                Radius for Lee-Sidford Ellipsoid (default is 0.5).
            g_lim : double, optional
                Constant for stopping gradient norm in gradient descent (default is 0.01).
            step_size : double, optional
                Constant for step size in gradient descent (default is 0.1).
            max_iter : int, optional
                Constant for maximum number of gradient descent iterations (default is 1000).
            err : double, optional
                Constant for error term in g^{-1}(x) (default is 1e-6).
            )doc",
            py::arg("r") = 0.5, py::arg("g_lim") = 0.01, 
            py::arg("step_size") = 0.1, py::arg("max_iter") = 1000, py::arg("err") = 1e-6);

}