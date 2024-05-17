#ifndef CONSDIKIN_HPP
#define CONSDIKIN_HPP

#include "SparseBarrierWalk.hpp"

class SparseDikinWalk : public SparseBarrierWalk{

    public:
        /**
         * @brief constructor for Dikin Walk class
         * @param r spread parameter
         * @param thin thin parameter
         * @param err error constant
         */
        SparseDikinWalk(double r, int thin = 1, double err = 1e-6) : SparseBarrierWalk(r, thin, err) {}

        /**
         * @brief generate weight (identity matrix)
         * @param x slack variable
         * @param A polytope constraint
         * @param k k values >= 0 constraint
         * @return SparseMatrixXd
         */
        SparseMatrixXd generateWeight(
            const VectorXd& x, 
            const SparseMatrixXd& A,
            int k
        ) override; 
    
    protected:

        /**
         * @brief Distribution constant
         * @param d polytope matrix 
         * @param n polytope vector
         * @return void
         */
        void setDistTerm(int d, int n) override;

};

#endif

