#ifndef SPARSE_LP_HPP
#define SPARSE_LP_HPP

#include "Common.hpp"
#include <ifopt/variable_set.h>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/problem.h>
#include <ifopt/ipopt_solver.h>
using namespace ifopt;

class SparseExVariables: public VariableSet{
    public:
        VectorXd x;

        SparseExVariables(int num_dim, string name, const VectorXd& init) : VariableSet(num_dim, name){
            x = init;
        } 

        void SetVariables(const VectorXd& val) override
        {
            x = val;
        }

        VectorXd GetValues() const override{
            return x;
        }

        VecBound GetBounds() const override{
            VecBound bounds(GetRows());
            for(int i = 0; i < x.rows(); i++){
                bounds.at(i) = NoBound;
            }
            return bounds;
        }

};

class SparseExConstraint1 : public ConstraintSet{
    public:
        SparseMatrixXd A;
        double err;
        SparseExConstraint1(int num_dim, string name, const SparseMatrixXd& A_param, double err_p = 1e-8) : ConstraintSet(num_dim, name){
            // A is d by n matrix
            A = A_param;
            err = err_p;
        }

        VectorXd GetValues() const override{
            VectorXd input = GetVariables()->GetComponent("var_set1")->GetValues();
            return A * input; 
        }

        VecBound GetBounds() const override{
            VecBound b(GetRows());
            for(int i = 0; i < A.rows(); i++){
                b.at(i) = Bounds(-inf, err);
            }
            return b;
        }
        
        void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{
            for(int k = 0; k < A.outerSize(); k++){
                for(SparseMatrixXd::InnerIterator it(A, k); it; ++it){
                    if (jac_block.coeffRef(it.row(), it.col()) != 0){
                        return; 
                    }
                    jac_block.coeffRef(it.row(), it.col()) = A.coeff(it.row(), it.col());
                }
            }
        }

};

class SparseExConstraint2 : public ConstraintSet{
    public:
        SparseMatrixXd A;
        VectorXd b;
        string name;
        double err;
   
        SparseExConstraint2(int num_dim, string name_, const SparseMatrixXd& A_param, const VectorXd& b_param, double err_p = 1e-8) : ConstraintSet(num_dim, name_){
            A = A_param;
            b = b_param;
            name = name_; 
            err = err_p;
        }

        VectorXd GetValues() const override{
            VectorXd input = GetVariables()->GetComponent(name)->GetValues();
            return A * input;
        }

        VecBound GetBounds() const override{
            VecBound bound(GetRows());
            for(int i = 0; i < b.rows(); i++){
                bound.at(i) = Bounds(b(i)-err,b(i)+err);
            }
            return bound;
        }

        void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{
            for(int k = 0; k < A.outerSize(); k++){
                for(SparseMatrixXd::InnerIterator it(A, k); it; ++it){
                    if (jac_block.coeffRef(it.row(), it.col()) != 0){
                        return; 
                    }
                    jac_block.coeffRef(it.row(), it.col()) = A.coeff(it.row(), it.col());
                }
            }

        }

};

class SparseExCost : public CostTerm{
    public:
        string name;

        SparseExCost(string name_) : CostTerm(name_) {
            name = name_; 
        }
        
        double GetCost() const override
        {
            VectorXd x = GetVariables()->GetComponent(name)->GetValues();
            return x(x.rows() - 1);
        }

        void FillJacobianBlock (string var_set, Jacobian& jac_block) const override{
            VectorXd x = GetVariables()->GetComponent(name)->GetValues();
            jac_block.coeffRef(0, x.rows() - 1) = 1.0; 

        }
    
};

#endif