// Copyright (C) 2018 Yixuan Qiu <yixuan.qiu@cos.name>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef ARNOLDI_OP_H
#define ARNOLDI_OP_H

#include <Eigen/Core>
#include <cmath>  // std::sqrt

#include <iostream>

namespace Spectra {


///
/// \ingroup Internals
/// @{
///

///
/// \defgroup Operators Operators
///
/// Different types of operators.
///

///
/// \ingroup Operators
///
/// Operators used in the Arnoldi factorization.
///
template <typename Scalar, typename OpType, typename BOpType>
class ArnoldiOp
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    OpType&  m_op;
    BOpType& m_Bop;
    Vector   m_cache;

public:
    ArnoldiOp(OpType* op, BOpType* Bop) :
        m_op(*op), m_Bop(*Bop), m_cache(op->rows())
    {}

    inline long int rows() const { return m_op.rows(); }

    // In generalized eigenvalue problem Ax=lambda*Bx, define the inner product to be <x, y> = x'By.
    // For regular eigenvalue problems, it is the usual inner product <x, y> = x'y

    // Compute <x, y> = x'By
    // x and y are two vectors
    template <typename Arg1, typename Arg2>
    Scalar inner_product(const Arg1& x, const Arg2& y)
    {
        m_Bop.mat_prod(y.data(), m_cache.data());
        return x.dot(m_cache);
    }

    // Compute res = <X, y> = X'By
    // X is a matrix, y is a vector, res is a vector
    template <typename Arg1, typename Arg2>
    void trans_product(const Arg1& x, const Arg2& y, Eigen::Ref<Vector> res)
    {
        m_Bop.mat_prod(y.data(), m_cache.data());
        res.noalias() = x.transpose() * m_cache;
    }

    // B-norm of a vector, ||x||_B = sqrt(x'Bx)
    template <typename Arg>
    Scalar norm(const Arg& x)
    {
        using std::sqrt;
        return sqrt(inner_product<Arg, Arg>(x, x));
    }

    // The "A" operator to generate the Krylov subspace
    inline void perform_op(const Scalar* x_in, Scalar* y_out)
    {
        m_op.perform_op(x_in, y_out);
    }
};



///
/// \ingroup Operators
///
/// Placeholder for the B-operator when \f$B = I\f$.
///
class IdentityBOp {};



///
/// \ingroup Operators
///
/// Partial specialization for the case \f$B = I\f$.
///
template <typename Scalar, typename OpType>
class ArnoldiOp<Scalar, OpType, IdentityBOp>
{
private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

    OpType& m_op;

public:
    ArnoldiOp<Scalar, OpType, IdentityBOp>(OpType* op, IdentityBOp* Bop) :
        m_op(*op)
    {}

    inline long int rows() const { return m_op.rows(); }

    // Compute <x, y> = x'y
    // x and y are two vectors
//    template <typename Arg1, typename Arg2>
//    Scalar inner_product(const Arg1& x, const Arg2& y) const
//    {
//        return x.dot(y);
//    }

    template <typename Arg1, typename Arg2>
    Scalar inner_product(const Arg1& x, const Arg2& y) const
    {
	Scalar sum = 0.;
        const int n_threads = omp_get_max_threads();
	std::vector<double> sum_thread(n_threads, 0.0);
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < x.size(); i++) {
            const int thread_id = omp_get_thread_num();
            sum_thread[thread_id] += x(i) * y(i);
        }
        for (int i = 0; i < n_threads; i++) sum += sum_thread[i];
        return sum;
    }

    // Compute res = <X, y> = X'y
    // X is a matrix, y is a vector, res is a vector
//    template <typename Arg1, typename Arg2>
//    void trans_product(const Arg1& x, const Arg2& y, Eigen::Ref<Vector> res) const
//    {
//        res.noalias() = x.transpose() * y;
//    }

    template <typename Arg1, typename Arg2>
    void trans_product(const Arg1& x, const Arg2& y, Eigen::Ref<Vector> res) const
    {
#pragma omp parallel for
        for (size_t i = 0; i < x.cols() ; i++) {
            //res(i) = inner_product<Arg2, Arg2>(x.col(i), y);
	    res(i) = x.col(i).dot(y);
	}

//	for (size_t i = 0; i < x.cols() ; i++) {
//	    res(i) = inner_product(x.col(i), y);
//	}
    }

    // B-norm of a vector. For regular eigenvalue problems it is simply the L2 norm
//    template <typename Arg>
//    Scalar norm(const Arg& x)
//    {
//        return x.norm();
//    }

    template <typename Arg>
    Scalar norm(const Arg& x)
    {
        return std::sqrt(inner_product(x, x));
    }

    // The "A" operator to generate the Krylov subspace
    inline void perform_op(const Scalar* x_in, Scalar* y_out)
    {
        m_op.perform_op(x_in, y_out);
    }

    // General matrix-vector multiplication in parallel
    // res = x * y. x is a matrix, y and z vectors
    template <typename Arg1, typename Arg2>
    void matvec_product(const Arg1& x, const Arg2& y, Eigen::Ref<Vector> res) const
    {
        res.setZero();
#pragma omp parallel for schedule(static)
        for (size_t k = 0; k < x.rows(); k++) {
            for (size_t j = 0; j < x.cols(); j++) {
                res(k) += x(k, j) * y(j);
            }
        }
    }
};

///
/// @}
///


} // namespace Spectra

#endif // ARNOLDI_OP_H
