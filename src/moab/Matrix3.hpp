/*
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 *
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 */

/**\file Matrix3.hpp
 *\author Jason Kraftcheck (kraftche@cae.wisc.edu)
 *\date 2006-07-18
 *\date 2012-08-2 Updated by rhl to be more generic. less code that does more!
 * TODO: Remove all 'inline' keywords as it is only a suggestion to the compiler
 * anyways, and it will ignore it or add it when it thinks its necessary.
 *\date 2016-08-03 Updated to use Eigen3 support underneath to improve performance
 */

#ifndef MOAB_MATRIX3_HPP
#define MOAB_MATRIX3_HPP

#include <iostream>
#include <iosfwd>
#include <limits>
#include <cmath>
#include <cassert>

#include "moab/MOABConfig.h"
#include "moab/ErrorHandler.hpp"
#include "moab/Util.hpp"
#include "moab/Types.hpp"
#include "moab/CartVect.hpp"

#ifdef MOAB_HAVE_EIGEN3

#if !defined( MOAB_HAVE_EIGEN3 ) && !defined( MOAB_HAVE_LAPACK )
// If we want unshifted QR iteration, uncomment below
// #define MOAB_EIGEN_DECOMPOSITION_UNSHIFTEDQR
#endif

#ifdef __GNUC__
// save diagnostic state
#pragma GCC diagnostic push
// turn off the specific warning. Can also use "-Wshadow"
#pragma GCC diagnostic ignored "-Wshadow"
#endif

#define EIGEN_DEFAULT_TO_ROW_MAJOR
#define EIGEN_INITIALIZE_MATRICES_BY_ZERO
// #define EIGEN_NO_STATIC_ASSERT
#include "Eigen/Dense"

#ifdef __GNUC__
// turn the warnings back on
#pragma GCC diagnostic pop
#endif

#else // ifdef MOAB_HAVE_LAPACK

#ifdef MOAB_HAVE_LAPACK

#if defined( MOAB_FC_FUNC_ )
#define MOAB_FC_WRAPPER MOAB_FC_FUNC_
#elif defined( MOAB_FC_FUNC )
#define MOAB_FC_WRAPPER MOAB_FC_FUNC
#else
#define MOAB_FC_WRAPPER( name, NAME ) name##_
#endif


// We will rely on LAPACK directly
#ifdef WIN32

// Should use second form below for windows but
// needed to do this to make it work.
// TODO: Need to clean this up
#define MOAB_dsyevd MOAB_FC_FUNC( dsyevd, DSYEVD )
#define MOAB_dgeev  MOAB_FC_FUNC( dgeev, DGEEV )

#else // ifndef WIN32

#define MOAB_dsyevd MOAB_FC_WRAPPER( dsyevd, DSYEVD )
#define MOAB_dgeev  MOAB_FC_WRAPPER( dgeev, DGEEV )

#endif // ifdef WIN32

extern "C" {

// Computes all eigenvalues and, optionally, eigenvectors of a
// real symmetric matrix A. If eigenvectors are desired, it uses a
// divide and conquer algorithm.
void MOAB_dsyevd( char* jobz,
                  char* uplo,
                  int* n,
                  double a[],
                  int* lda,
                  double w[],
                  double work[],
                  int* lwork,
                  int iwork[],
                  int* liwork,
                  int* info );


// Computes for an N-by-N real nonsymmetric matrix A, the
// eigenvalues and, optionally, the left and/or right eigenvectors.
void MOAB_dgeev( char* jobvl,
                 char* jobvr,
                 int* n,
                 double* a,
                 int* lda,
                 double* wr,
                 double* wi,
                 double* vl,
                 int* ldvl,
                 double* vr,
                 int* ldvr,
                 double* work,
                 int* lwork,
                 int* info );

}

#endif // ifdef MOAB_HAVE_LAPACK

#include <cstring>
#define MOAB_DMEMZERO( a, b ) memset( a, 0, ( b ) * sizeof( double ) )

#endif

namespace moab
{

namespace Matrix
{
    template < typename Matrix >
    inline Matrix mmult3( const Matrix& a, const Matrix& b )
    {
        return Matrix( a( 0, 0 ) * b( 0, 0 ) + a( 0, 1 ) * b( 1, 0 ) + a( 0, 2 ) * b( 2, 0 ),
                       a( 0, 0 ) * b( 0, 1 ) + a( 0, 1 ) * b( 1, 1 ) + a( 0, 2 ) * b( 2, 1 ),
                       a( 0, 0 ) * b( 0, 2 ) + a( 0, 1 ) * b( 1, 2 ) + a( 0, 2 ) * b( 2, 2 ),
                       a( 1, 0 ) * b( 0, 0 ) + a( 1, 1 ) * b( 1, 0 ) + a( 1, 2 ) * b( 2, 0 ),
                       a( 1, 0 ) * b( 0, 1 ) + a( 1, 1 ) * b( 1, 1 ) + a( 1, 2 ) * b( 2, 1 ),
                       a( 1, 0 ) * b( 0, 2 ) + a( 1, 1 ) * b( 1, 2 ) + a( 1, 2 ) * b( 2, 2 ),
                       a( 2, 0 ) * b( 0, 0 ) + a( 2, 1 ) * b( 1, 0 ) + a( 2, 2 ) * b( 2, 0 ),
                       a( 2, 0 ) * b( 0, 1 ) + a( 2, 1 ) * b( 1, 1 ) + a( 2, 2 ) * b( 2, 1 ),
                       a( 2, 0 ) * b( 0, 2 ) + a( 2, 1 ) * b( 1, 2 ) + a( 2, 2 ) * b( 2, 2 ) );
    }

    template < typename Matrix >
    inline const Matrix inverse( const Matrix& d )
    {
        const double det = 1.0 / determinant3( d );
        return inverse( d, det );
    }

    template < typename Vector, typename Matrix >
    inline Vector vector_matrix( const Vector& v, const Matrix& m )
    {
        return Vector( v[0] * m( 0, 0 ) + v[1] * m( 1, 0 ) + v[2] * m( 2, 0 ),
                       v[0] * m( 0, 1 ) + v[1] * m( 1, 1 ) + v[2] * m( 2, 1 ),
                       v[0] * m( 0, 2 ) + v[1] * m( 1, 2 ) + v[2] * m( 2, 2 ) );
    }

    template < typename Vector, typename Matrix >
    inline Vector matrix_vector( const Matrix& m, const Vector& v )
    {
        Vector res;
        res[0] = v[0] * m( 0, 0 ) + v[1] * m( 0, 1 ) + v[2] * m( 0, 2 );
        res[1] = v[0] * m( 1, 0 ) + v[1] * m( 1, 1 ) + v[2] * m( 1, 2 );
        res[2] = v[0] * m( 2, 0 ) + v[1] * m( 2, 1 ) + v[2] * m( 2, 2 );
        return res;
    }

}  // namespace Matrix

class Matrix3
{
  public:
    const static int size = 9;

  private:
#ifdef MOAB_HAVE_EIGEN3
    Eigen::Matrix3d _mat;
#else
    double _mat[size];
#endif

  public:
    // Default Constructor
    inline Matrix3()
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat.fill( 0.0 );
#else
        MOAB_DMEMZERO( _mat, Matrix3::size );
#endif
    }

#ifdef MOAB_HAVE_EIGEN3
    inline Matrix3( Eigen::Matrix3d mat ) : _mat( mat ) {}
#endif

    // TODO: Deprecate this.
    // Then we can go from three Constructors to one.
    inline Matrix3( double diagonal )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat << diagonal, 0.0, 0.0, 0.0, diagonal, 0.0, 0.0, 0.0, diagonal;
#else
        MOAB_DMEMZERO( _mat, Matrix3::size );
        _mat[0] = _mat[4] = _mat[8] = diagonal;
#endif
    }

    inline Matrix3( const CartVect& diagonal )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat << diagonal[0], 0.0, 0.0, 0.0, diagonal[1], 0.0, 0.0, 0.0, diagonal[2];
#else
        MOAB_DMEMZERO( _mat, Matrix3::size );
        _mat[0] = diagonal[0];
        _mat[4] = diagonal[1];
        _mat[8] = diagonal[2];
#endif
    }

    // TODO: not strictly correct as the Matrix3 object
    // is a double d[ 9] so the only valid model of T is
    // double, or any refinement (int, float)
    //*but* it doesn't really matter anything else
    // will fail to compile.
    inline Matrix3( const std::vector< double >& diagonal )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat << diagonal[0], 0.0, 0.0, 0.0, diagonal[1], 0.0, 0.0, 0.0, diagonal[2];
#else
        MOAB_DMEMZERO( _mat, Matrix3::size );
        _mat[0] = diagonal[0];
        _mat[4] = diagonal[1];
        _mat[8] = diagonal[2];
#endif
    }

    inline Matrix3( double v00,
                    double v01,
                    double v02,
                    double v10,
                    double v11,
                    double v12,
                    double v20,
                    double v21,
                    double v22 )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat << v00, v01, v02, v10, v11, v12, v20, v21, v22;
#else
        MOAB_DMEMZERO( _mat, Matrix3::size );
        _mat[0] = v00;
        _mat[1] = v01;
        _mat[2] = v02;
        _mat[3] = v10;
        _mat[4] = v11;
        _mat[5] = v12;
        _mat[6] = v20;
        _mat[7] = v21;
        _mat[8] = v22;
#endif
    }

    // Copy constructor
    Matrix3( const Matrix3& f )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat = f._mat;
#else
        memcpy( _mat, f._mat, size * sizeof( double ) );
#endif
    }

    // Weird constructors
    template < typename Vector >
    inline Matrix3( const Vector& row0, const Vector& row1, const Vector& row2, const bool isRow )
    {
#ifdef MOAB_HAVE_EIGEN3
        if( isRow )
        {
            _mat << row0[0], row0[1], row0[2], row1[0], row1[1], row1[2], row2[0], row2[1], row2[2];
        }
        else
        {
            _mat << row0[0], row1[0], row2[0], row0[1], row1[1], row2[1], row0[2], row1[2], row2[2];
        }
#else
        MOAB_DMEMZERO( _mat, Matrix3::size );
        if( isRow )
        {
            _mat[0] = row0[0];
            _mat[1] = row0[1];
            _mat[2] = row0[2];
            _mat[3] = row1[0];
            _mat[4] = row1[1];
            _mat[5] = row1[2];
            _mat[6] = row2[0];
            _mat[7] = row2[1];
            _mat[8] = row2[2];
        }
        else
        {
            _mat[0] = row0[0];
            _mat[1] = row1[0];
            _mat[2] = row2[0];
            _mat[3] = row0[1];
            _mat[4] = row1[1];
            _mat[5] = row2[1];
            _mat[6] = row0[2];
            _mat[7] = row1[2];
            _mat[8] = row2[2];
        }
#endif
    }

    /*
     * \deprecated { Use instead the constructor with explicit fourth argument, bool isRow, above }
     *
     */
    inline Matrix3( const double v[size] )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat << v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8];
#else
        memcpy( _mat, v, size * sizeof( double ) );
#endif
    }

    inline void copyto( double v[Matrix3::size] )
    {
#ifdef MOAB_HAVE_EIGEN3
        std::copy( _mat.data(), _mat.data() + size, v );
#else
        memcpy( v, _mat, size * sizeof( double ) );
#endif
    }

    inline Matrix3& operator=( const Matrix3& m )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat = m._mat;
#else
        memcpy( _mat, m._mat, size * sizeof( double ) );
#endif
        return *this;
    }

    inline Matrix3& operator=( const double v[size] )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat << v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8];
#else
        memcpy( _mat, v, size * sizeof( double ) );
#endif
        return *this;
    }

    inline double* operator[]( unsigned i )
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat.row( i ).data();
#else
        return &_mat[i * 3];  // Row Major
#endif
    }

    inline const double* operator[]( unsigned i ) const
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat.row( i ).data();
#else
        return &_mat[i * 3];
#endif
    }

    inline double& operator()( unsigned r, unsigned c )
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat( r, c );
#else
        return _mat[r * 3 + c];
#endif
    }

    inline double operator()( unsigned r, unsigned c ) const
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat( r, c );
#else
        return _mat[r * 3 + c];
#endif
    }

    inline double& operator()( unsigned i )
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat( i );
#else
        return _mat[i];
#endif
    }

    inline double operator()( unsigned i ) const
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat( i );
#else
        return _mat[i];
#endif
    }

    // get pointer to array of nine doubles
    inline double* array()
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat.data();
#else
        return _mat;
#endif
    }

    inline const double* array() const
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat.data();
#else
        return _mat;
#endif
    }

    inline Matrix3& operator+=( const Matrix3& m )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat += m._mat;
#else
        for( int i = 0; i < Matrix3::size; ++i )
            _mat[i] += m._mat[i];
#endif
        return *this;
    }

    inline Matrix3& operator-=( const Matrix3& m )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat -= m._mat;
#else
        for( int i = 0; i < Matrix3::size; ++i )
            _mat[i] -= m._mat[i];
#endif
        return *this;
    }

    inline Matrix3& operator*=( double s )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat *= s;
#else
        for( int i = 0; i < Matrix3::size; ++i )
            _mat[i] *= s;
#endif
        return *this;
    }

    inline Matrix3& operator/=( double s )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat /= s;
#else
        for( int i = 0; i < Matrix3::size; ++i )
            _mat[i] /= s;
#endif
        return *this;
    }

    inline Matrix3& operator*=( const Matrix3& m )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat *= m._mat;
#else
        // Uncomment below if you want point-wise multiplication instead (.*)
        // for (int i=0; i < Matrix3::size; ++i) _mat[i] *= m._mat[i];
        std::vector< double > dmat;
        dmat.assign( _mat, _mat + size );
        _mat[0] = dmat[0] * m._mat[0] + dmat[1] * m._mat[3] + dmat[2] * m._mat[6];
        _mat[1] = dmat[0] * m._mat[1] + dmat[1] * m._mat[4] + dmat[2] * m._mat[7];
        _mat[2] = dmat[0] * m._mat[2] + dmat[1] * m._mat[5] + dmat[2] * m._mat[8];
        _mat[3] = dmat[3] * m._mat[0] + dmat[4] * m._mat[3] + dmat[5] * m._mat[6];
        _mat[4] = dmat[3] * m._mat[1] + dmat[4] * m._mat[4] + dmat[5] * m._mat[7];
        _mat[5] = dmat[3] * m._mat[2] + dmat[4] * m._mat[5] + dmat[5] * m._mat[8];
        _mat[6] = dmat[6] * m._mat[0] + dmat[7] * m._mat[3] + dmat[8] * m._mat[6];
        _mat[7] = dmat[6] * m._mat[1] + dmat[7] * m._mat[4] + dmat[8] * m._mat[7];
        _mat[8] = dmat[6] * m._mat[2] + dmat[7] * m._mat[5] + dmat[8] * m._mat[8];
#endif
        return *this;
    }

    inline bool is_symmetric()
    {
        const double EPS = 1e-14;
#ifdef MOAB_HAVE_EIGEN3
        if( ( fabs( _mat( 1 ) - _mat( 3 ) ) < EPS ) && ( fabs( _mat( 2 ) - _mat( 6 ) ) < EPS ) &&
            ( fabs( _mat( 5 ) - _mat( 7 ) ) < EPS ) )
            return true;
#else
        if( ( fabs( _mat[1] - _mat[3] ) < EPS ) && ( fabs( _mat[2] - _mat[6] ) < EPS ) &&
            ( fabs( _mat[5] - _mat[7] ) < EPS ) )
            return true;
#endif
        else
            return false;
    }

    inline bool is_positive_definite()
    {
#ifdef MOAB_HAVE_EIGEN3
        double subdet6 = _mat( 1 ) * _mat( 5 ) - _mat( 2 ) * _mat( 4 );
        double subdet7 = _mat( 2 ) * _mat( 3 ) - _mat( 0 ) * _mat( 5 );
        double subdet8 = _mat( 0 ) * _mat( 4 ) - _mat( 1 ) * _mat( 3 );
        // Determinant:= d(6)*subdet6 + d(7)*subdet7 + d(8)*subdet8;
        const double det = _mat( 6 ) * subdet6 + _mat( 7 ) * subdet7 + _mat( 8 ) * subdet8;
        return _mat( 0 ) > 0 && subdet8 > 0 && det > 0;
#else
        double subdet6 = _mat[1] * _mat[5] - _mat[2] * _mat[4];
        double subdet7 = _mat[2] * _mat[3] - _mat[0] * _mat[5];
        double subdet8 = _mat[0] * _mat[4] - _mat[1] * _mat[3];
        // Determinant:= d(6)*subdet6 + d(7)*subdet7 + d(8)*subdet8;
        const double det = _mat[6] * subdet6 + _mat[7] * subdet7 + _mat[8] * subdet8;
        return _mat[0] > 0 && subdet8 > 0 && det > 0;
#endif
    }

#ifdef MOAB_HAVE_LAPACK
    template < typename Vector >
    inline ErrorCode eigen_decomposition_lapack( bool bisSymmetric, Vector& evals, Matrix3& evecs )
    {
        int info;
        /* Solve eigenproblem */
        double devreal[3], drevecs[9];
        if( !bisSymmetric )
        {
            double devimag[3], dlevecs[9], dwork[102];
            char dgeev_opts[2] = { 'N', 'V' };
            int N = 3, LWORK = 102, NL = 1, NR = N;
            std::vector< double > devmat(9);
            memcpy( devmat.data(), _mat, size * sizeof( double ) );
            // devmat.assign( _mat, _mat + size );
            MOAB_dgeev( &dgeev_opts[0], &dgeev_opts[1], &N, &devmat[0], &N, devreal, devimag, dlevecs, &NL, drevecs,
                        &NR, dwork, &LWORK, &info );
            // The result eigenvalues are ordered as high-->low
            evals[0]      = devreal[2];
            evals[1]      = devreal[1];
            evals[2]      = devreal[0];
            evecs._mat[0] = drevecs[6];
            evecs._mat[1] = drevecs[3];
            evecs._mat[2] = drevecs[0];
            evecs._mat[3] = drevecs[7];
            evecs._mat[4] = drevecs[4];
            evecs._mat[5] = drevecs[1];
            evecs._mat[6] = drevecs[8];
            evecs._mat[7] = drevecs[5];
            evecs._mat[8] = drevecs[2];
            std::cout << "DGEEV: Optimal work vector: dsize = " << dwork[0] << ".\n";
        }
        else
        {
            char dgeev_opts[2]      = { 'V', 'L' };
            const bool find_optimal = false;
            std::vector< int > iwork( 18 );
            std::vector< double > devmat( 9, 0.0 );
            std::vector< double > dwork( 38 );
            int N = 3, lwork = 38, liwork = 18;
            devmat[0] = _mat[0];
            devmat[1] = _mat[1];
            devmat[2] = _mat[2];
            devmat[4] = _mat[4];
            devmat[5] = _mat[5];
            devmat[8] = _mat[8];
            if( find_optimal )
            {
                int _lwork             = -1;
                int _liwork            = -1;
                double query_work_size = 0;
                int query_iwork_size   = 0;
                // Make an empty call to find the optimal work vector size
                MOAB_dsyevd( &dgeev_opts[0], &dgeev_opts[1], &N, NULL, &N, NULL, &query_work_size, &_lwork,
                             &query_iwork_size, &_liwork, &info );
                lwork = (int)query_work_size;
                dwork.resize( lwork );
                liwork = query_iwork_size;
                iwork.resize( liwork );
                std::cout << "DSYEVD: Optimal work vector: dsize = " << lwork << ", and isize = " << liwork << ".\n";
            }

            MOAB_dsyevd( &dgeev_opts[0], &dgeev_opts[1], &N, &devmat[0], &N, devreal, &dwork[0], &lwork, &iwork[0],
                         &liwork, &info );
            for( int i = 0; i < 9; ++i )
                drevecs[i] = devmat[i];
            // The result eigenvalues are ordered as low-->high, but vectors are in rows of A.
            evals[0]      = devreal[0];
            evals[1]      = devreal[1];
            evals[2]      = devreal[2];
            evecs._mat[0] = drevecs[0];
            evecs._mat[3] = drevecs[1];
            evecs._mat[6] = drevecs[2];
            evecs._mat[1] = drevecs[3];
            evecs._mat[4] = drevecs[4];
            evecs._mat[7] = drevecs[5];
            evecs._mat[2] = drevecs[6];
            evecs._mat[5] = drevecs[7];
            evecs._mat[8] = drevecs[8];
        }

        if( !info )
        {
            return MB_SUCCESS;
        }
        else
        {
            std::cout << "Failure in LAPACK_" << ( bisSymmetric ? "DSYEVD" : "DGEEV" )
                      << " call for eigen decomposition.\n";
            std::cout << "Failed with error = " << info << ".\n";
            return MB_FAILURE;
        }
    }
#endif

// #define TEST_NEWEIGEN
    template < typename Vector >
    inline ErrorCode eigen_decomposition( Vector& evals, Matrix3& evecs )
    {
#if defined( MOAB_HAVE_EIGEN3 )
        const bool bisSymmetric = this->is_symmetric();
        if( bisSymmetric )
        {
            Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d > eigensolver( this->_mat );
            if( eigensolver.info() != Eigen::Success ) return MB_FAILURE;
            const Eigen::SelfAdjointEigenSolver< Eigen::Matrix3d >::RealVectorType& e3evals = eigensolver.eigenvalues();
            evals[0]                                                                        = e3evals( 0 );
            evals[1]                                                                        = e3evals( 1 );
            evals[2]                                                                        = e3evals( 2 );
            evecs._mat = eigensolver.eigenvectors();  //.col(1)
            return MB_SUCCESS;
        }
        else
        {
            MB_CHK_SET_ERR( MB_FAILURE, "Unsymmetric matrix implementation with Eigen3 is currently not provided." );
            // Eigen::EigenSolver<Eigen::Matrix3d> eigensolver(this->_mat, true);
            // if (eigensolver.info() != Eigen::Success)
            //   return MB_FAILURE;
            // const Eigen::EigenSolver<Eigen::Matrix3d>::EigenvalueType& e3evals =
            // eigensolver.eigenvalues().real(); evals[0] = e3evals(0); evals[1] = e3evals(1);
            // evals[2] = e3evals(2); evecs._mat = eigensolver.eigenvectors().real(); //.col(1)
            // return MB_SUCCESS;
        }
#elif defined( MOAB_HAVE_LAPACK )
        const bool bisSymmetric = this->is_symmetric();
        return eigen_decomposition_lapack( bisSymmetric, evals, evecs );
        // return eigen_decomposition_native( evals, evecs );
        // eigen_decomposition_lapack( bisSymmetric, evals, evecs );
        // Vector evals_native;
        // Matrix3 evecs_native;
        // eigen_decomposition_native( evals_native, evecs_native );
        // // std::cout << "LAPACK: " << evals << " NATIVE: " << evals_native << " Error: " << evals-evals_native << std::endl;
        // evals = evals_native;
        // evecs = evecs_native;
        // return MB_SUCCESS;
#elif defined( TEST_NEWEIGEN )
        // Helper function to compute the absolute maximum off-diagonal element
        auto maxOffDiagonal = [] (double* matrix) -> std::pair<int, int> {
            // int p = 0, q = 0;
            std::pair<int,int> pq;
            double maxVal = 0.0;
            for (int i = 0; i < 3; ++i) {
                for (int j = i + 1; j < 3; ++j) {
                    double absVal = std::abs(matrix[i * 3 + j]);
                    if (absVal > maxVal) {
                        maxVal = absVal;
                        pq     = std::pair< int, int >( i, j );
                    }
                }
            }
            return pq;
        };

        // Jacobi rotation for 3x3 matrix
        auto jacobiRotate = [](double* matrix, double* eigenvectors, int p, int q) {
            double app = matrix[p * 3 + p];
            double aqq = matrix[q * 3 + q];
            double apq = matrix[p * 3 + q];

            // Compute the Jacobi rotation
            double phi = 0.5 * std::atan2(2 * apq, aqq - app);
            double c = std::cos(phi);
            double s = std::sin(phi);

            // Perform rotation
            for (int i = 0; i < 3; ++i) {
                double aip = matrix[i * 3 + p];
                double aiq = matrix[i * 3 + q];
                matrix[i * 3 + p] = c * aip - s * aiq;
                matrix[i * 3 + q] = s * aip + c * aiq;
            }
            for (int i = 0; i < 3; ++i) {
                double aip = matrix[p * 3 + i];
                double aiq = matrix[q * 3 + i];
                matrix[p * 3 + i] = c * aip - s * aiq;
                matrix[q * 3 + i] = s * aip + c * aiq;
            }

            // Update eigenvectors
            for (int i = 0; i < 3; ++i) {
                double vip = eigenvectors[i * 3 + p];
                double viq = eigenvectors[i * 3 + q];
                eigenvectors[i * 3 + p] = c * vip - s * viq;
                eigenvectors[i * 3 + q] = s * vip + c * viq;
            }
        };

        // Initialize the matrix and eigenvector matrix
    evecs = Matrix3::identity();
    Matrix3 A = *this; // make a copy

    const double tolerance = std::numeric_limits<double>::epsilon();
    const int maxIterations = 100;
    bool converged = false;

    for (int iter = 0; iter < maxIterations; ++iter) {
        // Find the largest off-diagonal element
        auto pqpair = maxOffDiagonal( A._mat );
        if( std::abs( A._mat[pqpair.first * 3 + pqpair.second] ) < tolerance )
        {
            converged = true;
            break;  // Converged
        }

        // Apply Jacobi rotation
        jacobiRotate( A._mat, evecs._mat, pqpair.first, pqpair.second );
    }

    // Eigenvalues are now on the diagonal
    evals[0] = A._mat[0];
    evals[1] = A._mat[4];
    evals[2] = A._mat[8];

    // compute the transpose of the vector
    // evecs.transpose_inplace();

    Vector eval_native;
    Matrix3 evec_native;
    eigen_decomposition_native(eval_native, evec_native);

    if( !converged )
    {
        std::cout << "Eigen decomposition not converged \n";
        A.print( std::cout );
        std::cout << evals << std::endl;
        evecs.print( std::cout );

        std::cout << eval_native << std::endl;
        evec_native.print( std::cout );
        evals = eval_native;
        evecs = evec_native;
        return MB_FAILURE;
    }

    return (converged ? MB_SUCCESS : MB_FAILURE);

#elif defined(MOAB_EIGEN_DECOMPOSITION_UNSHIFTEDQR)

    // QR decomposition using Gram-Schmidt process
    auto qrDecomposition = []( const Matrix3& A, Matrix3& Q, Matrix3& R ) {
        // Matrix3 q, a;
        moab::CartVect q[3], a[3];

        // Extract columns of A as vectors and copy to a
        for( int i = 0; i < 3; ++i )
        {
            a[i] = A.vcol<moab::CartVect>( i );
            // a[i][0] = A._mat[0 + i];
            // a[i][1] = A._mat[3 + i];
            // a[i][2] = A._mat[6 + i];
        }

        // Gram-Schmidt process
        for( int i = 0; i < 3; ++i )
        {
            q[i] = a[i];
            for( int j = 0; j < i; ++j )
            {
                double dot = q[j] % a[i];  // Dot product
                q[i] -= dot * q[j];
            }
            q[i].normalize();
        }

        // Form Q matrix (from q vectors)
        for( int i = 0; i < 3; ++i )
        {
            Q._mat[0 + i] = q[i][0];
            Q._mat[3 + i] = q[i][1];
            Q._mat[6 + i] = q[i][2];
        }

        // Form R matrix
        for( int i = 0; i < 3; ++i )
        {
            for( int j = i; j < 3; ++j )
            {
                R( i, j ) = q[i] % a[j];  // Dot product
            }
        }

    };

    // QR algorithm to compute eigenvalues and eigenvectors
    Matrix3& A = *this;
    // qrAlgorithm
    int maxIter = 100;
    double tol  = 1e-8;
    bool converged = false;
    {
        Matrix3 Q( 0.0 );  // Orthogonal matrix
        Matrix3 R( 0.0 );  // Upper triangular matrix

        // Initialize eigenvectors as identity matrix
        evecs = Matrix3( 1.0, 0.0, 0.0,  // Identity matrix
                         0.0, 1.0, 0.0,  // Orthogonal
                         0.0, 0.0, 1.0 );

        for( int iter = 0; iter < maxIter; ++iter )
        {
            // QR decomposition of A
            qrDecomposition( A, Q, R );

            // Compute A = R * Q
            A = R * Q;

            // Update eigenvectors (accumulate Q)
            evecs = evecs * Q;

            // Check convergence by looking at the off-diagonal elements
            double offDiagonal = std::abs( A( 1, 0 ) ) + std::abs( A( 2, 0 ) ) + std::abs( A( 2, 1 ) );
            if( offDiagonal < tol )
            {
                converged = true;
                break;
            }
        }

        // Extract eigenvalues (diagonal elements of A)
        evals[0] = A( 0, 0 );
        evals[1] = A( 1, 1 );
        evals[2] = A( 2, 2 );
    }

    if( !converged )
    {
        std::cout << "Eigen decomposition not converged \n";
        A.print(std::cout);
        std::cout << evals << std::endl;
        evecs.print(std::cout);
        exit( 1 );
    }

    return ( converged ? MB_SUCCESS : MB_FAILURE );
#else
    return eigen_decomposition_native( evals, evecs );
#endif
    }

    // Function to perform a Jacobi rotation
    void jacobiRotate( moab::Matrix3& A, moab::Matrix3& V, int p, int q )
    {
        if( A(p,q) == 0 ) return;

        double theta, t, c, s;
        theta = ( A(q,q) - A(p,p) ) / ( 2.0 * A(p,q) );
        t     = ( theta >= 0 ) ? 1.0 / ( theta + std::sqrt( 1.0 + theta * theta ) )
                               : 1.0 / ( theta - std::sqrt( 1.0 + theta * theta ) );
        c     = 1.0 / std::sqrt( 1.0 + t * t );
        s     = t * c;

        double app = A(p,p), aqq = A(q,q), apq = A(p,q);
        A(p,p) = c * c * app - 2.0 * c * s * apq + s * s * aqq;
        A(q,q) = s * s * app + 2.0 * c * s * apq + c * c * aqq;
        A(p,q) = A(q,p) = 0.0;  // Zero the off-diagonal element

        // Update other elements
        for( int i = 0; i < 3; i++ )
        {
            if( i != p && i != q )
            {
                double aip = A(i,p), aiq = A(i,q);
                A(i,p) = A(p,i) = c * aip - s * aiq;
                A(i,q) = A(q,i) = s * aip + c * aiq;
            }
        }

        // Update the eigenvector matrix
        for( int i = 0; i < 3; i++ )
        {
            double vip = V(i,p), viq = V(i,q);
            V(i,p) = c * vip - s * viq;
            V(i,q) = s * vip + c * viq;
        }
    }

    // Jacobi's method to find eigenvalues and eigenvectors of a symmetric 3x3 matrix
    moab::ErrorCode jacobiEigenDecomposition( moab::Matrix3& A, moab::CartVect& eigenvalues, moab::Matrix3& eigenvectors )
    {
        constexpr double EPSILON = 1e-14;  // Tolerance for stopping the iteration
        constexpr int maxiters   = 500;
        constexpr int n = 3;

        // Initialize the eigenvector matrix as the identity matrix
        eigenvectors = moab::Matrix3::identity();

        // Iterate to apply Jacobi rotations
        bool converged = false;
        for( int iter = 0; iter < maxiters; ++iter )
        {
            // Find the largest off-diagonal element
            int p = 0, q = 1;
            double maxOffDiag = std::abs( A(p,q) );
            for( int i = 0; i < n; ++i )
            {
                for( int j = i + 1; j < n; ++j )
                {
                    if( std::abs( A(i,j) ) > maxOffDiag )
                    {
                        maxOffDiag = std::abs( A(i,j) );
                        p          = i;
                        q          = j;
                    }
                }
            }

            // If the largest off-diagonal element is smaller than tolerance, stop
            if( maxOffDiag < EPSILON )
            {
                converged = true;
                break;
            }

            // Apply Jacobi rotation to zero out A[p][q]
            jacobiRotate( A, eigenvectors, p, q );
        }

        // The diagonal elements of A are the eigenvalues
        for( int i = 0; i < n; ++i )
        {
            // Extract eigenvalues (diagonal elements of A)
            eigenvalues[i] = A( i, i );
            eigenvectors.col(i).normalize();
        }

        if( !converged )
        {
            std::cerr << "Jacobi method did not converge within " << maxiters << " iterations\n";
        }
        else
        {
            auto sortIndices = []( const CartVect& vec ) -> std::array< int, 3 > {
                // Initialize the index array with values 0, 1, 2
                std::array< int, 3 > indices = { 0, 1, 2 };

                // Sort the indices based on the values in the original vector
                std::sort( indices.begin(), indices.end(), [&]( int i, int j ) { return vec[i] < vec[j]; } );

                return indices;
            };

            auto newIndices = sortIndices( eigenvalues );

            // eigenvectors.transpose_inplace();
            CartVect teigenvalues = eigenvalues;
            Matrix3 teigenvectors = eigenvectors;
            for( int i = 0; i < 3; i++ )
            {
                eigenvalues[i]       = teigenvalues[newIndices[i]];
                eigenvectors( i, 0 ) = teigenvectors( i, newIndices[0] );
                eigenvectors( i, 1 ) = teigenvectors( i, newIndices[1] );
                eigenvectors( i, 2 ) = teigenvectors( i, newIndices[2] );
            }
        }

        return ( converged ? MB_SUCCESS : MB_FAILURE );
    }

    template < typename Vector >
    inline ErrorCode eigen_decomposition_native( Vector& evals, Matrix3& evecs )
    {
        // return eigen_decomposition_native_arnoldiqr( evals, evecs );
        // return eigen_decomposition_native_shiftedqr( evals, evecs );
        // eigen_decomposition_native_analytical( evals, evecs );

        Matrix3 A = *this;
        return jacobiEigenDecomposition( A, evals, evecs );

        // if( !converged )
        // if (false)
        // {
        //     std::cout.precision( 16 );

        //     std::cout << "\nMatrix \n";
        //     this->print( std::cout );

        //     Vector evals_lap;
        //     Matrix3 evecs_lap;
        //     eigen_decomposition_lapack( true, evals_lap, evecs_lap );
        //     std::cout << "\nEigenvalues and eigenvectors LAPACK \n";
        //     std::cout << evals_lap << std::endl;
        //     std::cout << "\nEigenvectors \n";
        //     evecs_lap.print( std::cout );

        //     std::cout << "\nEigenvalues and eigenvectors Native \n";
        //     std::cout << evals << std::endl;
        //     std::cout << "\nEigenvectors \n";
        //     evecs.print( std::cout );
        //     // A.print( std::cout );
        // }
        return MB_SUCCESS;
    }

    template < typename Vector >
    inline ErrorCode eigen_decomposition_native_arnoldiqr( Vector& evals, Matrix3& evecs )
    {
        // Function to perform the Modified Gram-Schmidt orthogonalization
        // auto modifiedGramSchmidt = []( const std::vector< moab::CartVect >& Q, moab::CartVect& v, int k ) {
        //     for( int i = 0; i < k; ++i )
        //     {
        //         double dot_product = Q[i] % v;  // Dot product between v and Q[i]
        //         v -= dot_product * Q[i];        // Orthogonalize
        //     }
        // };

        // Arnoldi Iteration to approximate the eigenvalues and eigenvectors of a 3x3 matrix
        auto arnoldiIteration = []( const moab::Matrix3& A, std::vector< moab::CartVect >& Q, moab::Matrix3& H ) {
            constexpr int maxIter = 3;

            // Initialize the first vector (random start, here using a simple vector)
            moab::CartVect v0( 1.0, 1.0, 1.0 );
            v0.normalize();
            Q.push_back( v0 );  // Orthonormal basis, first vector

            // Iterate and build the Hessenberg matrix H
            for( int k = 0; k < maxIter; ++k )
            {
                // Multiply the matrix A by the current vector Q[k]
                // moab::CartVect v = A * Q[k];
                auto v = Matrix::matrix_vector( A, Q[k] );

                // Modified Gram-Schmidt orthogonalization
                // modifiedGramSchmidt( Q, v, k + 1 );
                for( int im = 0; im < k+1; ++im )
                {
                    double dot_product = Q[im] % v;  // Dot product between v and Q[i]
                    v -= dot_product * Q[im];        // Orthogonalize
                }

                // Compute the norm of v
                double norm_v = v.length();

                // Store in Hessenberg matrix
                if( k < 2 )
                {
                    H( k + 1, k ) = norm_v;
                }

                if( norm_v > 1e-9 )
                {
                    v /= norm_v;       // Normalize the new vector
                    Q.push_back( v );  // Add it to the orthonormal basis
                }

                // Fill the upper part of Hessenberg matrix H
                for( int i = 0; i <= k; ++i )
                    H( i, k ) = Q[i] % Matrix::matrix_vector( A, Q[k] );
            }
        };

        // QR decomposition using Gram-Schmidt process
        auto qrDecomposition = []( const Matrix3& A, Matrix3& Q, Matrix3& R ) {
            moab::CartVect q[3], a[3];

            // Extract columns of A as vectors and copy to a
            for( int i = 0; i < 3; ++i )
            {
                a[i] = A.vcol< moab::CartVect >( i );
            }

            // Gram-Schmidt process
            for( int i = 0; i < 3; ++i )
            {
                q[i] = a[i];
                for( int j = 0; j < i; ++j )
                {
                    double dot = q[j] % a[i];  // Dot product
                    q[i] -= dot * q[j];
                }
                q[i].normalize();
            }

            // Form Q matrix (from q vectors)
            for( int i = 0; i < 3; ++i )
            {
                Q( 0, i ) = q[i][0];
                Q( 1, i ) = q[i][1];
                Q( 2, i ) = q[i][2];
            }

            // Form R matrix
            for( int i = 0; i < 3; ++i )
                for( int j = i; j < 3; ++j )
                    R( i, j ) = q[i] % a[j];  // Dot product
        };

        auto qrAlgorithm =
            [&qrDecomposition]( const moab::Matrix3& matrix, moab::CartVect& eigenvalues,
                                moab::Matrix3& eigenvectors ) {
                // Initialize eigenvectors as identity matrix
                eigenvectors = moab::Matrix3( 1.0 );

                // QR algorithm to compute eigenvalues and eigenvectors
                // Matrix3 A = *this;
                // Start the shifted-QR iteration
                int maxIter    = 5000;
                double tol     = 1e-8;
                bool converged = false;

                moab::Matrix3 Q( 0.0 );  // Orthogonal matrix
                moab::Matrix3 R( 0.0 );  // Upper triangular matrix

                double d, mu;
                moab::Matrix3 A = matrix, Ap = matrix;
                for( int iter = 0; iter < maxIter; ++iter )
                {
                    // Wilkinson shift: use bottom-right 2x2 submatrix to compute shift
                    d = ( A( 1, 1 ) - A( 2, 2 ) ) / 2.0;
                    if( fabs( d ) < tol )  // possible multiplicity; skip shift
                        mu = 0.0;          // disable shift
                    else
                        mu = A( 2, 2 ) - ( d / std::abs( d ) ) * A( 2, 1 ) * A( 2, 1 ) /
                                             ( std::abs( d ) + std::sqrt( d * d + A( 2, 1 ) * A( 2, 1 ) ) );


                    //  mu = 0.0;

                    // Apply the shift
                    moab::Matrix3 A_shifted = A - Matrix3( mu );

                    // QR decomposition of A
                    qrDecomposition( A_shifted, Q, R );

                    // Compute A = R * Q + mu * I (undo the shift)
                    A = R * Q + Matrix3( mu );

                    // Update eigenvectors (accumulate Q)
                    eigenvectors = eigenvectors * Q;

                    // Check convergence by looking at the off-diagonal elements
                    Ap -= A;
                    double offDiagonal =
                        std::sqrt( Ap( 1, 0 ) * Ap( 1, 0 ) + Ap( 2, 0 ) * Ap( 2, 0 ) + Ap( 2, 1 ) * Ap( 2, 1 ) );
                    if( offDiagonal < tol )
                    {
                        converged = true;
                        break;
                    }
                    else
                        Ap = A;
                }

                auto sortIndices = []( const Vector& vec ) -> std::array< int, 3 > {
                    // Initialize the index array with values 0, 1, 2
                    std::array< int, 3 > indices = { 0, 1, 2 };

                    // Sort the indices based on the values in the original vector
                    std::sort( indices.begin(), indices.end(), [&]( int i, int j ) { return vec[i] < vec[j]; } );

                    return indices;
                };

                // Extract eigenvalues (diagonal elements of A)
                eigenvalues[0] = A( 0, 0 );
                eigenvalues[1] = A( 1, 1 );
                eigenvalues[2] = A( 2, 2 );

                if( converged )
                {
                    auto newIndices = sortIndices( eigenvalues );

                    // eigenvectors.transpose_inplace();
                    Vector teigenvalues   = eigenvalues;
                    Matrix3 teigenvectors = eigenvectors;
                    for( int i = 0; i < 3; i++ )
                    {
                        eigenvalues[i]       = teigenvalues[newIndices[i]];
                        eigenvectors( i, 0 ) = teigenvectors( i, newIndices[0] );
                        eigenvectors( i, 1 ) = teigenvectors( i, newIndices[1] );
                        eigenvectors( i, 2 ) = teigenvectors( i, newIndices[2] );
                    }
                }

                return converged;
            };

        auto computeEigenvectors = []( const std::vector< moab::CartVect >& Q, const moab::Matrix3& eigenvectors_H,
                                       moab::Matrix3& eigenvectors ) {
            // Eigenvectors of the original matrix A are Q * y (where y is an eigenvector of H)
            for( int i = 0; i < 3; ++i )
            {
                moab::CartVect eigenvector( 0.0, 0.0, 0.0 );
                for( int j = 0; j < 3; ++j )
                {
                    eigenvector += eigenvectors_H( j, i ) * Q[j];  // Map back to original space
                }
                eigenvector.normalize();
                eigenvectors( 0, i ) = eigenvector[0];
                eigenvectors( 1, i ) = eigenvector[1];
                eigenvectors( 2, i ) = eigenvector[2];
            }
        };

        // Vectors to hold the orthonormal basis (Q) and Hessenberg matrix (H)
        std::vector< moab::CartVect > Q;
        moab::Matrix3 H( 0.0 );  // Initialize the Hessenberg matrix to zero

        // Perform Arnoldi iteration
        arnoldiIteration( *this, Q, H );

        // Compute eigenvalues from Hessenberg matrix
        moab::Matrix3 eigenvectors_h( 1.0 );
        bool converged = qrAlgorithm( H, evals, eigenvectors_h );

        if( !converged )
        {
            std::cout << "\nEigen decomposition not converged \n";
            std::cout.precision( 16 );

            std::cout << "\nMatrix \n";
            this->print( std::cout );

            Vector evals_lap;
            Matrix3 evecs_lap;
            eigen_decomposition_lapack( true, evals_lap, evecs_lap );
            std::cout << "\nEigenvalues and eigenvectors LAPACK \n";
            std::cout << evals_lap << std::endl;
            std::cout << "\nEigenvectors \n";
            evecs_lap.print( std::cout );

            std::cout << "\nEigenvalues and eigenvectors Native \n";
            std::cout << evals << std::endl;
            std::cout << "\nEigenvectors \n";
            evecs.print( std::cout );
            // A.print( std::cout );
            return MB_FAILURE;
        }
        else
        {
            // std::cout << "\nEigenvalues of the matrix:" << std::endl;
            // std::cout << eigenvalues[0] << " " << eigenvalues[1] << " " << eigenvalues[2] << std::endl;

            // Compute the eigenvectors of the original matrix
            computeEigenvectors( Q, eigenvectors_h, evecs );
            return MB_SUCCESS;
        }
    }
//
    template < typename Vector >
    inline ErrorCode eigen_decomposition_native_analytical( Vector& evals, Matrix3& evecs )
    {
        // taken from Eigen3: struct direct_selfadjoint_eigenvalues<SolverType,3,false>
        typedef double Scalar;
        Matrix3 m = *this;

        // Shift the matrix to the mean eigenvalue and map the matrix coefficients to [-1:1] to avoid over- and underflow.
        Scalar shift = m.trace() / Scalar( 3 );
        // TODO Avoid this copy. Currently it is necessary to suppress bogus values when determining maxCoeff and for computing the eigenvectors later
        Matrix3 scaledMat = m;
        // scaledMat.diagonal().array() -= shift;
        scaledMat._mat[0] -= shift;
        scaledMat._mat[4] -= shift;
        scaledMat._mat[8] -= shift;
        Scalar scale = std::numeric_limits<double>::min();
        // scale = scaledMat.cwiseAbs().maxCoeff();
        for( unsigned i = 0; i < 9; ++i )
            if( scale < fabs( scaledMat._mat[i] ) ) scale = fabs( scaledMat._mat[i] );
        if( scale > 0 ) scaledMat /= scale;  // TODO for scale==0 we could save the remaining operations

        // Now let us analytically compute the eigenvalues
        const Scalar s_inv3  = Scalar( 1 ) / Scalar( 3 );
        const Scalar s_sqrt3 = std::sqrt( Scalar( 3 ) );

        // The characteristic equation is x^3 - c2*x^2 + c1*x - c0 = 0.  The
        // eigenvalues are the roots to this equation, all guaranteed to be
        // real-valued, because the matrix is symmetric.
        Scalar c0 = m( 0, 0 ) * m( 1, 1 ) * m( 2, 2 ) + Scalar( 2 ) * m( 1, 0 ) * m( 2, 0 ) * m( 2, 1 ) -
                    m( 0, 0 ) * m( 2, 1 ) * m( 2, 1 ) - m( 1, 1 ) * m( 2, 0 ) * m( 2, 0 ) -
                    m( 2, 2 ) * m( 1, 0 ) * m( 1, 0 );
        Scalar c1 = m( 0, 0 ) * m( 1, 1 ) - m( 1, 0 ) * m( 1, 0 ) + m( 0, 0 ) * m( 2, 2 ) - m( 2, 0 ) * m( 2, 0 ) +
                    m( 1, 1 ) * m( 2, 2 ) - m( 2, 1 ) * m( 2, 1 );
        Scalar c2 = m( 0, 0 ) + m( 1, 1 ) + m( 2, 2 );

#define MATRIX_MAXI(x,y) (x < y ? y : x)
        // Construct the parameters used in classifying the roots of the equation
        // and in solving the equation for the roots in closed form.
        Scalar c2_over_3 = c2 * s_inv3;
        Scalar a_over_3  = ( c2 * c2_over_3 - c1 ) * s_inv3;
        a_over_3         = MATRIX_MAXI( a_over_3, Scalar( 0 ) );

        Scalar half_b = Scalar( 0.5 ) * ( c0 + c2_over_3 * ( Scalar( 2 ) * c2_over_3 * c2_over_3 - c1 ) );

        Scalar q = a_over_3 * a_over_3 * a_over_3 - half_b * half_b;
        q        = MATRIX_MAXI( q, Scalar( 0 ) );

        // Compute the eigenvalues by solving for the roots of the polynomial.
        Scalar rho = sqrt( a_over_3 );
        Scalar theta =
            atan2( sqrt( q ), half_b ) * s_inv3;  // since sqrt(q) > 0, atan2 is in [0, pi] and theta is in [0, pi/3]
        Scalar cos_theta = cos( theta );
        Scalar sin_theta = sin( theta );
        // roots are already sorted, since cos is monotonically decreasing on [0, pi]
        evals[0] = c2_over_3 - rho * ( cos_theta + s_sqrt3 * sin_theta );  // == 2*rho*cos(theta+2pi/3)
        evals[1] = c2_over_3 - rho * ( cos_theta - s_sqrt3 * sin_theta );  // == 2*rho*cos(theta+ pi/3)
        evals[2] = c2_over_3 + Scalar( 2 ) * rho * cos_theta;

        // Rescale back to the original size.
        evals *= scale;
        evals += shift;
        return moab::MB_SUCCESS;

        /// Now let us compute eigenvectors

        // auto extract_kernel = []( moab::Matrix3& mat, CartVect& res, CartVect& representative ) -> bool {
        //     int i0;
        //     // Find non-zero column i0 (by construction, there must exist a non zero coefficient on the diagonal):
        //     mat.diagonal().cwiseAbs().maxCoeff( &i0 );
        //     // mat.col(i0) is a good candidate for an orthogonal vector to the current eigenvector,
        //     // so let's save it:
        //     representative = mat.col( i0 );
        //     Scalar n0, n1;
        //     moab::CartVect c0, c1;
        //     n0 = ( c0 = representative.cross( mat.col( ( i0 + 1 ) % 3 ) ) ).squaredNorm();
        //     n1 = ( c1 = representative.cross( mat.col( ( i0 + 2 ) % 3 ) ) ).squaredNorm();
        //     if( n0 > n1 )
        //         res = c0 / std::sqrt( n0 );
        //     else
        //         res = c1 / std::sqrt( n1 );

        //     return true;
        // };

        // Function to extract an orthogonal vector based on the matrix and representative vector
        // bool extract_kernel( moab::Matrix3 & mat, moab::CartVect & res, moab::CartVect & representative )
        auto extract_kernel = []( moab::Matrix3& mat, CartVect& res, CartVect& representative ) -> bool {
            int i0 = -1;

            // Find the index of the largest absolute diagonal element of the matrix
            moab::CartVect diag( mat( 0, 0 ), mat( 1, 1 ), mat( 2, 2 ) );  // Get the diagonal elements
            double maxVal = std::abs( diag[0] );
            i0            = 0;

            for( int i = 1; i < 3; ++i )
            {
                if( std::abs( diag[i] ) > maxVal )
                {
                    maxVal = std::abs( diag[i] );
                    i0     = i;
                }
            }

            // Extract the column of the matrix corresponding to the index `i0`
            representative = mat.vcol< moab::CartVect >( i0 );

            // Variables to store cross product results and their norms
            moab::CartVect c0, c1;
            double n0, n1;

            // Compute the cross products between the representative vector and other matrix columns
            moab::CartVect col1 = mat.vcol< moab::CartVect >( ( i0 + 1 ) % 3 );
            moab::CartVect col2 = mat.vcol< moab::CartVect >( ( i0 + 2 ) % 3 );

            // Cross product of representative with col1 and col2
            c0 = representative * col1;
            c1 = representative * col2;

            // Compute the squared norms of the cross products
            n0 = c0.length_squared();
            n1 = c1.length_squared();

            // Select the larger norm cross product and normalize the result
            if( n0 > n1 )
            {
                res = c0 / std::sqrt( n0 );
            }
            else
            {
                res = c1 / std::sqrt( n1 );
            }

            return true;
        };

        if( ( fabs( evals[2] - evals[0] ) ) <= std::numeric_limits< double >::min() )
        {
            // All three eigenvalues are numerically the same
            evecs = moab::Matrix3( 1.0 );
            }
            else
            {
                Matrix3 tmp;
                tmp = scaledMat;

                // Compute the eigenvector of the most distinct eigenvalue
                Scalar d0 = evals[2] - evals[1];
                Scalar d1 = evals[1] - evals[0];
                int k( 0 ), l( 2 );
                if( d0 > d1 )
                {
                    int t = l;
                    l     = k;
                    k     = t;
                    d0    = d1;
                }

                // Compute the eigenvector of index k
                {
                    // tmp.diagonal().array() -= evals( k );
                    tmp._mat[0] -= evals[k];
                    tmp._mat[4] -= evals[k];
                    tmp._mat[8] -= evals[k];
                    // By construction, 'tmp' is of rank 2, and its kernel corresponds to the respective eigenvector.
                    moab::CartVect col_k = evecs.col( k );
                    moab::CartVect col_l = evecs.col( l );
                    extract_kernel( tmp, col_k, col_l );
                    evecs( 0, l ) = col_l[0];
                    evecs( 1, l ) = col_l[1];
                    evecs( 2, l ) = col_l[2];
                    evecs( 0, k ) = col_k[0];
                    evecs( 1, k ) = col_k[1];
                    evecs( 2, k ) = col_k[2];
                }

                // Compute eigenvector of index l
                if( d0 <= 2 * std::numeric_limits< double >::min() * d1 )
                {
                    // If d0 is too small, then the two other eigenvalues are numerically the same,
                    // and thus we only have to ortho-normalize the near orthogonal vector we saved above.
                    evecs.col( l ) -= ( evecs.col( k ) % evecs.col( l ) ) * evecs.col( l );
                    evecs.col( l ).normalize();
                }
                else
                {
                    tmp = scaledMat;
                    // tmp.diagonal().array() -= evals( l );
                    tmp._mat[0] -= evals[l];
                    tmp._mat[4] -= evals[l];
                    tmp._mat[8] -= evals[l];

                    moab::CartVect dummy;
                    moab::CartVect col_l = evecs.col( l );
                    extract_kernel( tmp, col_l, dummy );
                    evecs( 0, l ) = col_l[0];
                    evecs( 1, l ) = col_l[1];
                    evecs( 2, l ) = col_l[2];
                }
                // Compute last eigenvector from the other two
                CartVect lev = ( evecs.col( 2 ) * evecs.col( 0 ) );
                lev.normalize(); evecs.col( 1 ) = lev;
            }

            // Rescale back to the original size.
            evals *= scale;
            evals += shift;

            return moab::MB_SUCCESS;
        }

        template < typename Vector >
        inline ErrorCode eigen_decomposition_native_shiftedqr( Vector & evals, Matrix3 & evecs )
        {
            // QR decomposition using Gram-Schmidt process
            auto qrDecomposition = []( const Matrix3& A, Matrix3& Q, Matrix3& R ) {
                moab::CartVect q[3], a[3];

                // Extract columns of A as vectors and copy to a
                for( int i = 0; i < 3; ++i )
                {
                    a[i] = A.vcol< moab::CartVect >( i );
                    // a[i][0] = A._mat[0 + i];
                    // a[i][1] = A._mat[3 + i];
                    // a[i][2] = A._mat[6 + i];
                }

                // Gram-Schmidt process
                for( int i = 0; i < 3; ++i )
                {
                    q[i] = a[i];
                    for( int j = 0; j < i; ++j )
                    {
                        double dot = q[j] % a[i];  // Dot product
                        q[i] -= dot * q[j];
                    }
                    q[i].normalize();
                }

                // Form Q matrix (from q vectors)
                for( int i = 0; i < 3; ++i )
                {
                    Q._mat[0 + i] = q[i][0];
                    Q._mat[3 + i] = q[i][1];
                    Q._mat[6 + i] = q[i][2];
                }

                // Form R matrix
                for( int i = 0; i < 3; ++i )
                    for( int j = i; j < 3; ++j )
                        R( i, j ) = q[i] % a[j];  // Dot product
            };

            const Matrix3 I = moab::Matrix3( 1.0 );

            // QR algorithm to compute eigenvalues and eigenvectors
            Matrix3 A = *this;
            // Start the shifted-QR iteration
            int maxIter    = 2000;
            double tol     = 1e-10;
            bool converged = false;
            {
                Matrix3 Q( 0.0 );  // Orthogonal matrix
                Matrix3 R( 0.0 );  // Upper triangular matrix

                // Initialize eigenvectors as identity matrix
                evecs = I;

                double d, mu;
                Matrix3 Ap = A;
                for( int iter = 0; iter < maxIter; ++iter )
                {
                    // Wilkinson shift: use bottom-right 2x2 submatrix to compute shift
                    d = ( A( 1, 1 ) - A( 2, 2 ) ) / 2.0;
                    if( fabs( d ) < tol )  // possible multiplicity; skip shift
                        mu = 0.0;          // disable shift
                    else
                        mu = A( 2, 2 ) - ( d / std::abs( d ) ) * A( 2, 1 ) * A( 2, 1 ) /
                                             ( std::abs( d ) + std::sqrt( d * d + A( 2, 1 ) * A( 2, 1 ) ) );

                    //  mu = 0.0;

                    // Apply the shift
                    moab::Matrix3 A_shifted = A - Matrix3( mu );

                    // QR decomposition of A
                    qrDecomposition( A_shifted, Q, R );

                    // Compute A = R * Q + mu * I (undo the shift)
                    A = R * Q + Matrix3( mu );

                    // Update eigenvectors (accumulate Q)
                    evecs = evecs * Q;

                    // Check convergence by looking at the off-diagonal elements
                    Ap -= A;
                    double offDiagonal =
                        std::sqrt( Ap( 1, 0 ) * Ap( 1, 0 ) + Ap( 2, 0 ) * Ap( 2, 0 ) + Ap( 2, 1 ) * Ap( 2, 1 ) );
                    if( offDiagonal < tol )
                    {
                        converged = true;
                        break;
                    }
                    else
                        Ap = A;
                }

                auto sortIndices = []( const Vector& vec ) -> std::array< int, 3 > {
                    // Initialize the index array with values 0, 1, 2
                    std::array< int, 3 > indices = { 0, 1, 2 };

                    // Sort the indices based on the values in the original vector
                    std::sort( indices.begin(), indices.end(), [&]( int i, int j ) { return vec[i] < vec[j]; } );

                    return indices;
                };

                // Extract eigenvalues (diagonal elements of A)
                evals[0] = A( 0, 0 );
                evals[1] = A( 1, 1 );
                evals[2] = A( 2, 2 );

                auto newIndices = sortIndices( evals );

                // evecs.transpose_inplace();

                Vector tevals  = evals;
                Matrix3 tevecs = evecs;
                for( int i = 0; i < 3; i++ )
                {
                    evals[i]      = tevals[newIndices[i]];
                    evecs( i, 0 ) = tevecs( i, newIndices[0] );
                    evecs( i, 1 ) = tevecs( i, newIndices[1] );
                    evecs( i, 2 ) = tevecs( i, newIndices[2] );

                    // evecs( 0, i ) = tevecs( newIndices[0], i );
                    // evecs( 1, i ) = tevecs( newIndices[1], i );
                    // evecs( 2, i ) = tevecs( newIndices[2], i );
                }
            }

            if( !converged )
            {
                std::cout << "\nEigen decomposition not converged \n";
                std::cout.precision( 16 );

                std::cout << "\nMatrix \n";
                this->print( std::cout );

                Vector evals_lap;
                Matrix3 evecs_lap;
                eigen_decomposition_lapack( true, evals_lap, evecs_lap );
                std::cout << "\nEigenvalues and eigenvectors LAPACK \n";
                std::cout << evals_lap << std::endl;
                std::cout << "\nEigenvectors \n";
                evecs_lap.print( std::cout );

                std::cout << "\nEigenvalues and eigenvectors Native \n";
                std::cout << evals << std::endl;
                std::cout << "\nEigenvectors \n";
                evecs.print( std::cout );
                // A.print( std::cout );
                return MB_FAILURE;
            }

            return ( converged ? MB_SUCCESS : MB_FAILURE );
    }

    inline static Matrix3 identity()
    {
        // Identity matrix
        return Matrix3( 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 );
    }

    inline void transpose_inplace()
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat.transposeInPlace();
#else
        Matrix3 mtmp( *this );
        _mat[1] = mtmp._mat[3];
        _mat[3] = mtmp._mat[1];
        _mat[2] = mtmp._mat[6];
        _mat[6] = mtmp._mat[2];
        _mat[5] = mtmp._mat[7];
        _mat[7] = mtmp._mat[5];
#endif
    }

    inline Matrix3 transpose() const
    {
#ifdef MOAB_HAVE_EIGEN3
        return Matrix3( _mat.transpose() );
#else
        Matrix3 mtmp( *this );
        mtmp._mat[1] = _mat[3];
        mtmp._mat[3] = _mat[1];
        mtmp._mat[2] = _mat[6];
        mtmp._mat[6] = _mat[2];
        mtmp._mat[5] = _mat[7];
        mtmp._mat[7] = _mat[5];
        return mtmp;
#endif
    }

    template < typename Vector >
    inline void copycol( int index, Vector& vol )
    {
#ifdef MOAB_HAVE_EIGEN3
        _mat.col( index ).swap( vol );
#else
        switch( index )
        {
            case 0:
                _mat[0] = vol[0];
                _mat[3] = vol[1];
                _mat[6] = vol[2];
                break;
            case 1:
                _mat[1] = vol[0];
                _mat[4] = vol[1];
                _mat[7] = vol[2];
                break;
            case 2:
                _mat[2] = vol[0];
                _mat[5] = vol[1];
                _mat[8] = vol[2];
                break;
        }
#endif
    }

    inline void swapcol( int srcindex, int destindex )
    {
        assert( srcindex < Matrix3::size );
        assert( destindex < Matrix3::size );
#ifdef MOAB_HAVE_EIGEN3
        _mat.col( srcindex ).swap( _mat.col( destindex ) );
#else
        CartVect svol = this->vcol< CartVect >( srcindex );
        CartVect dvol = this->vcol< CartVect >( destindex );
        switch( srcindex )
        {
            case 0:
                _mat[0] = dvol[0];
                _mat[3] = dvol[1];
                _mat[6] = dvol[2];
                break;
            case 1:
                _mat[1] = dvol[0];
                _mat[4] = dvol[1];
                _mat[7] = dvol[2];
                break;
            case 2:
                _mat[2] = dvol[0];
                _mat[5] = dvol[1];
                _mat[8] = dvol[2];
                break;
        }
        switch( destindex )
        {
            case 0:
                _mat[0] = svol[0];
                _mat[3] = svol[1];
                _mat[6] = svol[2];
                break;
            case 1:
                _mat[1] = svol[0];
                _mat[4] = svol[1];
                _mat[7] = svol[2];
                break;
            case 2:
                _mat[2] = svol[0];
                _mat[5] = svol[1];
                _mat[8] = svol[2];
                break;
        }
#endif
    }

    template < typename Vector >
    inline Vector vcol( int index ) const
    {
        assert( index < Matrix3::size );
#ifdef MOAB_HAVE_EIGEN3
        return _mat.col( index );
#else
        switch( index )
        {
            case 0:
                return Vector( _mat[0], _mat[3], _mat[6] );
            case 1:
                return Vector( _mat[1], _mat[4], _mat[7] );
            case 2:
                return Vector( _mat[2], _mat[5], _mat[8] );
        }
        return Vector( 0.0 );
#endif
    }

    inline void colscale( int index, double scale )
    {
        assert( index < Matrix3::size );
#ifdef MOAB_HAVE_EIGEN3
        _mat.col( index ) *= scale;
#else
        switch( index )
        {
            case 0:
                _mat[0] *= scale;
                _mat[3] *= scale;
                _mat[6] *= scale;
                break;
            case 1:
                _mat[1] *= scale;
                _mat[4] *= scale;
                _mat[7] *= scale;
                break;
            case 2:
                _mat[2] *= scale;
                _mat[5] *= scale;
                _mat[8] *= scale;
                break;
        }
#endif
    }

    inline void rowscale( int index, double scale )
    {
        assert( index < Matrix3::size );
#ifdef MOAB_HAVE_EIGEN3
        _mat.row( index ) *= scale;
#else
        switch( index )
        {
            case 0:
                _mat[0] *= scale;
                _mat[1] *= scale;
                _mat[2] *= scale;
                break;
            case 1:
                _mat[3] *= scale;
                _mat[4] *= scale;
                _mat[5] *= scale;
                break;
            case 2:
                _mat[6] *= scale;
                _mat[7] *= scale;
                _mat[8] *= scale;
                break;
        }
#endif
    }

    inline CartVect col( int index ) const
    {
        assert( index < Matrix3::size );
#ifdef MOAB_HAVE_EIGEN3
        Eigen::Vector3d mvec = _mat.col( index );
        return CartVect( mvec[0], mvec[1], mvec[2] );
#else
        switch( index )
        {
            case 0:
                return CartVect( _mat[0], _mat[3], _mat[6] );
            case 1:
                return CartVect( _mat[1], _mat[4], _mat[7] );
            case 2:
                return CartVect( _mat[2], _mat[5], _mat[8] );
        }
        return CartVect( 0.0 );
#endif
    }

    inline CartVect row( int index ) const
    {
        assert( index < Matrix3::size );
#ifdef MOAB_HAVE_EIGEN3
        Eigen::Vector3d mvec = _mat.row( index );
        return CartVect( mvec[0], mvec[1], mvec[2] );
#else
        switch( index )
        {
            case 0:
                return CartVect( _mat[0], _mat[1], _mat[2] );
            case 1:
                return CartVect( _mat[3], _mat[4], _mat[5] );
            case 2:
                return CartVect( _mat[6], _mat[7], _mat[8] );
        }
        return CartVect( 0.0 );
#endif
    }

    inline CartVect diagonal( ) const
    {
#ifdef MOAB_HAVE_EIGEN3
        return CartVect( _mat( 0, 0 ), _mat( 1, 1 ), _mat( 2, 2 ) );
#else
        return CartVect( _mat[0], _mat[4], _mat[8] );
#endif
    }

    inline double trace() const
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat.trace();
#else
        return _mat[0] + _mat[4] + _mat[8];
#endif
    }

    friend Matrix3 operator+( const Matrix3& a, const Matrix3& b );
    friend Matrix3 operator-( const Matrix3& a, const Matrix3& b );
    friend Matrix3 operator*( const Matrix3& a, const Matrix3& b );

    inline double determinant() const
    {
#ifdef MOAB_HAVE_EIGEN3
        return _mat.determinant();
#else
        return ( _mat[0] * _mat[4] * _mat[8] + _mat[1] * _mat[5] * _mat[6] + _mat[2] * _mat[3] * _mat[7] -
                 _mat[0] * _mat[5] * _mat[7] - _mat[1] * _mat[3] * _mat[8] - _mat[2] * _mat[4] * _mat[6] );
#endif
    }

    inline Matrix3 inverse() const
    {
#ifdef MOAB_HAVE_EIGEN3
        return Matrix3( _mat.inverse() );
#else
        // return Matrix::compute_inverse( *this, this->determinant() );
        Matrix3 m( 0.0 );
        const double d_determinant = 1.0 / this->determinant();
        m._mat[0]                  = d_determinant * ( _mat[4] * _mat[8] - _mat[5] * _mat[7] );
        m._mat[1]                  = d_determinant * ( _mat[2] * _mat[7] - _mat[8] * _mat[1] );
        m._mat[2]                  = d_determinant * ( _mat[1] * _mat[5] - _mat[4] * _mat[2] );
        m._mat[3]                  = d_determinant * ( _mat[5] * _mat[6] - _mat[8] * _mat[3] );
        m._mat[4]                  = d_determinant * ( _mat[0] * _mat[8] - _mat[6] * _mat[2] );
        m._mat[5]                  = d_determinant * ( _mat[2] * _mat[3] - _mat[5] * _mat[0] );
        m._mat[6]                  = d_determinant * ( _mat[3] * _mat[7] - _mat[6] * _mat[4] );
        m._mat[7]                  = d_determinant * ( _mat[1] * _mat[6] - _mat[7] * _mat[0] );
        m._mat[8]                  = d_determinant * ( _mat[0] * _mat[4] - _mat[3] * _mat[1] );
        return m;
#endif
    }

    inline bool invert()
    {
        bool invertible = false;
        double d_determinant;
#ifdef MOAB_HAVE_EIGEN3
        Eigen::Matrix3d invMat;
        _mat.computeInverseAndDetWithCheck( invMat, d_determinant, invertible );
        if( !Util::is_finite( d_determinant ) ) return false;
        _mat = invMat;
        return invertible;
#else
        d_determinant = this->determinant();
        if( d_determinant > 1e-13 ) invertible = true;
        d_determinant = 1.0 / d_determinant;  // invert the determinant
        std::vector< double > _m;
        _m.assign( _mat, _mat + size );
        _mat[0]      = d_determinant * ( _m[4] * _m[8] - _m[5] * _m[7] );
        _mat[1]      = d_determinant * ( _m[2] * _m[7] - _m[8] * _m[1] );
        _mat[2]      = d_determinant * ( _m[1] * _m[5] - _m[4] * _m[2] );
        _mat[3]      = d_determinant * ( _m[5] * _m[6] - _m[8] * _m[3] );
        _mat[4]      = d_determinant * ( _m[0] * _m[8] - _m[6] * _m[2] );
        _mat[5]      = d_determinant * ( _m[2] * _m[3] - _m[5] * _m[0] );
        _mat[6]      = d_determinant * ( _m[3] * _m[7] - _m[6] * _m[4] );
        _mat[7]      = d_determinant * ( _m[1] * _m[6] - _m[7] * _m[0] );
        _mat[8]      = d_determinant * ( _m[0] * _m[4] - _m[3] * _m[1] );
#endif
        return invertible;
    }

    // Calculate determinant of 2x2 submatrix composed of the
    // elements not in the passed row or column.
    inline double subdet( int r, int c ) const
    {
        assert( r >= 0 && c >= 0 );
        if( r < 0 || c < 0 ) return DBL_MAX;
#ifdef MOAB_HAVE_EIGEN3
        const int r1 = ( r + 1 ) % 3, r2 = ( r + 2 ) % 3;
        const int c1 = ( c + 1 ) % 3, c2 = ( c + 2 ) % 3;
        return _mat( r1, c1 ) * _mat( r2, c2 ) - _mat( r1, c2 ) * _mat( r2, c1 );
#else
        const int r1 = Matrix3::size * ( ( r + 1 ) % 3 ), r2 = Matrix3::size * ( ( r + 2 ) % 3 );
        const int c1 = ( c + 1 ) % 3, c2 = ( c + 2 ) % 3;
        return _mat[r1 + c1] * _mat[r2 + c2] - _mat[r1 + c2] * _mat[r2 + c1];
#endif
    }

    inline void print( std::ostream& s ) const
    {
#ifdef MOAB_HAVE_EIGEN3
        s << "| " << _mat( 0 ) << " " << _mat( 1 ) << " " << _mat( 2 ) << " | " << _mat( 3 ) << " " << _mat( 4 ) << " "
          << _mat( 5 ) << " | " << _mat( 6 ) << " " << _mat( 7 ) << " " << _mat( 8 ) << " |";
#else
        s << "| " << _mat[0] << " " << _mat[1] << " " << _mat[2] << " | " << _mat[3] << " " << _mat[4] << " " << _mat[5]
          << " | " << _mat[6] << " " << _mat[7] << " " << _mat[8] << " |";
#endif
    }

};  // class Matrix3

template < typename Vector >
inline Matrix3 outer_product( const Vector& u, const Vector& v )
{
    return Matrix3( u[0] * v[0], u[0] * v[1], u[0] * v[2], u[1] * v[0], u[1] * v[1], u[1] * v[2], u[2] * v[0],
                    u[2] * v[1], u[2] * v[2] );
}

inline Matrix3 operator+( const Matrix3& a, const Matrix3& b )
{
#ifdef MOAB_HAVE_EIGEN3
    return Matrix3( a._mat + b._mat );
#else
    Matrix3 s( a );
    for( int i = 0; i < Matrix3::size; ++i )
        s( i ) += b._mat[i];
    return s;
#endif
}

inline Matrix3 operator-( const Matrix3& a, const Matrix3& b )
{
#ifdef MOAB_HAVE_EIGEN3
    return Matrix3( a._mat - b._mat );
#else
    Matrix3 s( a );
    for( int i = 0; i < Matrix3::size; ++i )
        s( i ) -= b._mat[i];
    return s;
#endif
}

inline Matrix3 operator*( const Matrix3& a, const Matrix3& b )
{
#ifdef MOAB_HAVE_EIGEN3
    return Matrix3( a._mat * b._mat );
#else
    return Matrix::mmult3( a, b );
#endif
}

template < typename T >
inline std::vector< T > operator*( const Matrix3& m, const std::vector< T >& v )
{
    return moab::Matrix::matrix_vector( m, v );
}

template < typename T >
inline std::vector< T > operator*( const std::vector< T >& v, const Matrix3& m )
{
    return moab::Matrix::vector_matrix( v, m );
}

inline CartVect operator*( const Matrix3& m, const CartVect& v )
{
    return moab::Matrix::matrix_vector( m, v );
}

inline CartVect operator*( const CartVect& v, const Matrix3& m )
{
    return moab::Matrix::vector_matrix( v, m );
}

}  // namespace moab

#ifdef MOAB_HAVE_LAPACK
#undef MOAB_DMEMZERO
#endif

#ifndef MOAB_MATRIX3_OPERATORLESS
#define MOAB_MATRIX3_OPERATORLESS
inline std::ostream& operator<<( std::ostream& s, const moab::Matrix3& m )
{
    return s << "| " << m( 0, 0 ) << " " << m( 0, 1 ) << " " << m( 0, 2 ) << " | " << m( 1, 0 ) << " " << m( 1, 1 )
             << " " << m( 1, 2 ) << " | " << m( 2, 0 ) << " " << m( 2, 1 ) << " " << m( 2, 2 ) << " |";
}
#endif  // MOAB_MATRIX3_OPERATORLESS

#endif  // MOAB_MATRIX3_HPP
