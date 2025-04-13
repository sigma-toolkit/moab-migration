#ifndef MBA_MBA_HPP
#define MBA_MBA_HPP

/*
The MIT License

Copyright (c) 2015 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   mba/mba.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Multilevel B-spline interpolation.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <array>
#include <map>
#include <list>
#include <memory>
#include <algorithm>
#include <numeric>
#include <functional>
#include <type_traits>
#include <utility>
#include <iterator>
#include <stdexcept>

#include <cmath>
#include <cassert>

#include "moab/Remapping/TempestRemapper.hpp"

namespace mba
{

using std::ptrdiff_t;
using std::size_t;

template < int N >
using point = std::array< double, N >;

template < int N >
using index = std::array< size_t, N >;

namespace detail
{

    const int32_t initial_tag = 42000;
    const int32_t dense_tag   = 42001;
    const int32_t sparse_tag  = 42002;

    template < class Cond, class Msg >
    void precondition( const Cond& cond, const Msg& msg )
    {
        if( !static_cast< bool >( cond ) ) throw std::runtime_error( msg );
    }

    // Compile-time N^M
    template < size_t N, size_t M >
    struct power : std::integral_constant< size_t, N * power< N, M - 1 >::value >
    {
    };

    template < size_t N >
    struct power< N, 0 > : std::integral_constant< size_t, 1 >
    {
    };

    // N-dimensional dense matrix
    template < class T, int N >
    class multi_array
    {
        static_assert( N > 0, "Wrong number of dimensions" );

      public:
        multi_array() {}

        multi_array( index< N > n )
        {
            init( n );
        }

        void resize( index< N > n )
        {
            init( n );
        }

        size_t size() const
        {
            return buf.size();
        }

        T operator()( index< N > i ) const
        {
            return buf[idx( i )];
        }

        T& operator()( index< N > i )
        {
            return buf[idx( i )];
        }

        T operator[]( size_t i ) const
        {
            return buf[i];
        }

        T& operator[]( size_t i )
        {
            return buf[i];
        }

        const T* data() const
        {
            return buf.data();
        }

        T* data()
        {
            return buf.data();
        }

        void read( std::ifstream& f )
        {
            precondition( static_cast< bool >( f.read( (char*)sizes.data(), sizeof( sizes ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)stride.data(), sizeof( stride ) ) ), "file i/o error" );

            size_t n = 1;
            for( size_t i = 0; i < N; ++i )
                n *= sizes[i];

            buf.resize( n );
            precondition( static_cast< bool >( f.read( (char*)buf.data(), sizeof( T ) * n ) ), "file i/o error" );
        }

        void write( std::ofstream& f ) const
        {
            precondition( static_cast< bool >( f.write( (char*)sizes.data(), sizeof( sizes ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)stride.data(), sizeof( stride ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)buf.data(), sizeof( T ) * buf.size() ) ),
                          "file i/o error" );
        }

      private:
        std::array< int, N > sizes;
        std::array< int, N > stride;
        std::vector< T > buf;

        void init( index< N > n )
        {
            size_t s = 1;

            for( int d = N - 1; d >= 0; --d )
            {
                sizes[d]  = n[d];
                stride[d] = s;
                s *= n[d];
            }

            buf.resize( s );
        }

        size_t idx( index< N > i ) const
        {
            size_t p = 0;
            for( int d = 0; d < N; ++d )
            {
#ifndef NDEBUG
                if( i[d] >= static_cast<size_t>(sizes[d]) )
                    throw std::out_of_range( "Index is out of range in mba::multi_array" );
#endif
                p += stride[d] * i[d];
            }
            return p;
        }
    };

    /// N-dimensional grid iterator (nested loop with variable depth).
    template < unsigned NDim >
    class grid_iterator
    {
      public:
        explicit grid_iterator( const index< NDim >& dims ) : N( dims ), i{ 0 }, done( i == N ), idx( 0 ) {}

        explicit grid_iterator( size_t dim ) : i{ 0 }, done( dim == 0 ), idx( 0 )
        {
            for( auto& v : N )
                v = dim;
        }

        size_t operator[]( size_t d ) const
        {
            return i[d];
        }

        const index< NDim >& operator*() const
        {
            return i;
        }

        size_t position() const
        {
            return idx;
        }

        grid_iterator& operator++()
        {
            done = true;
            for( size_t d = NDim; d--; )
            {
                if( ++i[d] < N[d] )
                {
                    done = false;
                    break;
                }
                i[d] = 0;
            }

            ++idx;

            return *this;
        }

        operator bool() const
        {
            return !done;
        }

      private:
        index< NDim > N, i;
        bool done;
        size_t idx;
    };

    template < class I, class F >
    struct transform_iterator
    {
        typedef typename I::iterator_category iterator_category;
        typedef ptrdiff_t difference_type;
        typedef size_t size_type;
        typedef decltype( std::declval< F >()( *std::declval< I >() ) ) value_type;
        typedef value_type* pointer;
        typedef value_type& reference;

        I i;
        F f;

        transform_iterator( I i, F f ) : i( i ), f( f ) {}
        transform_iterator( const transform_iterator& other ) : i( other.i ), f( other.f ) {}

        auto operator*() const -> decltype( f( *i ) )
        {
            return f( *i );
        }

        transform_iterator& operator++()
        {
            ++i;
            return *this;
        }

        bool operator!=( const transform_iterator& other ) const
        {
            return i != other.i;
        }

        ptrdiff_t operator-( const transform_iterator& other ) const
        {
            return i - other.i;
        }
    };

    template < class I, class F >
    transform_iterator< I, F > make_transform_iterator( I i, F f )
    {
        return transform_iterator< I, F >( i, f );
    }

    template < typename T, size_t N >
    std::array< T, N > operator+( std::array< T, N > a, const std::array< T, N >& b )
    {
        for( size_t i = 0; i < N; ++i )
            a[i] += b[i];
        return a;
    }

    template < typename T, size_t N, typename C >
    std::array< T, N > operator-( std::array< T, N > a, C b )
    {
        for( auto& v : a )
            v -= b;
        return a;
    }

    template < typename T, size_t N, typename C >
    std::array< T, N > operator*( std::array< T, N > a, C b )
    {
        for( auto& v : a )
            v *= b;
        return a;
    }

    // Value of k-th B-Spline basic function at t.
    inline double Bspline( size_t k, double t )
    {
        assert( 0 <= t && t < 1 );
        assert( k < 4 );

        switch( k )
        {
            case 0:
                return ( t * ( t * ( -t + 3 ) - 3 ) + 1 ) / 6;
            case 1:
                return ( t * t * ( 3 * t - 6 ) + 4 ) / 6;
            case 2:
                return ( t * ( t * ( -3 * t + 3 ) + 3 ) + 1 ) / 6;
            case 3:
                return t * t * t / 6;
            default:
                return 0;
        }
    }

    // Checks if p is between lo and hi
    template < typename T, size_t N >
    bool boxed( const std::array< T, N >& lo, const std::array< T, N >& p, const std::array< T, N >& hi )
    {
        for( unsigned i = 0; i < N; ++i )
        {
            if( p[i] < lo[i] || p[i] > hi[i] ) return false;
        }
        return true;
    }

    inline double safe_divide( double a, double b )
    {
        return b == 0.0 ? 0.0 : a / b;
    }

    template < unsigned NDim >
    class control_lattice
    {
      public:
        virtual ~control_lattice() {}

        virtual double operator()( const point< NDim >& p ) const = 0;

        virtual void report( std::ostream& ) const = 0;
        virtual int32_t tag() const                = 0;
        virtual void write( std::ofstream& ) const = 0;

        template < class CooIter, class ValIter >
        double residual( CooIter coo_begin, CooIter coo_end, ValIter val_begin ) const
        {
            double res = 0.0;

            CooIter p = coo_begin;
            ValIter v = val_begin;

            for( ; p != coo_end; ++p, ++v )
            {
                ( *v ) -= ( *this )( *p );
                res = std::max( res, std::abs( *v ) );
            }

            return res;
        }
    };

    template < unsigned NDim >
    class initial_approximation : public control_lattice< NDim >
    {
      public:
        initial_approximation( std::function< double( const point< NDim >& ) > f ) : f( f ) {}

        double operator()( const point< NDim >& p ) const
        {
            return f( p );
        }

        void report( std::ostream& os ) const
        {
            os << "initial approximation";
        }

        int32_t tag() const
        {
            return detail::initial_tag;
        }

        void write( std::ofstream& ) const {}

      private:
        std::function< double( const point< NDim >& ) > f;
    };

    template < unsigned NDim >
    class control_lattice_dense : public control_lattice< NDim >
    {
      public:
        template < class CooIter, class ValIter >
        control_lattice_dense( const point< NDim >& coo_min,
                               const point< NDim >& coo_max,
                               index< NDim > grid_size,
                               CooIter coo_begin,
                               CooIter coo_end,
                               ValIter val_begin )
            : cmin( coo_min ), cmax( coo_max ), grid( grid_size )
        {
            static_assert( std::is_same< typename std::iterator_traits< CooIter >::iterator_category,
                                         std::random_access_iterator_tag >::value,
                           "CooIter should be a random access iterator" );
            static_assert( std::is_same< typename std::iterator_traits< ValIter >::iterator_category,
                                         std::random_access_iterator_tag >::value,
                           "ValIter should be a random access iterator" );

            for( unsigned i = 0; i < NDim; ++i )
            {
                hinv[i] = ( grid[i] - 1 ) / ( cmax[i] - cmin[i] );
                cmin[i] -= 1 / hinv[i];
                grid[i] += 2;
            }

            multi_array< double, NDim > delta( grid );
            multi_array< double, NDim > omega( grid );

            std::fill( delta.data(), delta.data() + delta.size(), 0.0 );
            std::fill( omega.data(), omega.data() + omega.size(), 0.0 );

            phi.resize( grid );
            std::fill( phi.data(), phi.data() + phi.size(), 0.0 );

            ptrdiff_t n = std::distance( coo_begin, coo_end );
            size_t m    = phi.size();

#pragma omp parallel
            {
                multi_array< double, NDim > t_delta( grid );
                multi_array< double, NDim > t_omega( grid );

                std::fill( t_delta.data(), t_delta.data() + t_delta.size(), 0.0 );
                std::fill( t_omega.data(), t_omega.data() + t_omega.size(), 0.0 );

#pragma omp for schedule( dynamic, 1 )
                for( ptrdiff_t l = 0; l < n; ++l )
                {
                    auto p = coo_begin[l];
                    auto v = val_begin[l];

                    if( !boxed( coo_min, p, coo_max ) ) continue;

                    index< NDim > i;
                    point< NDim > s;

                    for( unsigned d = 0; d < NDim; ++d )
                    {
                        double u = ( p[d] - cmin[d] ) * hinv[d];
                        i[d]     = floor( u ) - 1;
                        s[d]     = u - floor( u );
                    }

                    std::array< double, power< 4, NDim >::value > w;
                    double sum_w2 = 0.0;

                    for( grid_iterator< NDim > d( 4 ); d; ++d )
                    {
                        double prod = 1.0;
                        for( unsigned k = 0; k < NDim; ++k )
                            prod *= Bspline( d[k], s[k] );

                        w[d.position()] = prod;
                        sum_w2 += prod * prod;
                    }

                    for( grid_iterator< NDim > d( 4 ); d; ++d )
                    {
                        double w1  = w[d.position()];
                        double w2  = w1 * w1;
                        double phi = v * w1 / sum_w2;

                        index< NDim > j = i + ( *d );

                        t_delta( j ) += w2 * phi;
                        t_omega( j ) += w2;
                    }
                }

                {
                    for( size_t i = 0; i < m; ++i )
                    {
#pragma omp atomic
                        delta[i] += t_delta[i];
#pragma omp atomic
                        omega[i] += t_omega[i];
                    }
                }
            }

            for( size_t i = 0; i < m; ++i )
            {
                phi[i] = safe_divide( delta[i], omega[i] );
            }
        }

        control_lattice_dense( std::ifstream& f )
        {
            precondition( static_cast< bool >( f.read( (char*)cmin.data(), sizeof( cmin ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)cmax.data(), sizeof( cmax ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)hinv.data(), sizeof( hinv ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)grid.data(), sizeof( grid ) ) ), "file i/o error" );

            phi.read( f );
        }

        void write( std::ofstream& f ) const
        {
            precondition( static_cast< bool >( f.write( (char*)cmin.data(), sizeof( cmin ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)cmax.data(), sizeof( cmax ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)hinv.data(), sizeof( hinv ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)grid.data(), sizeof( grid ) ) ), "file i/o error" );

            phi.write( f );
        }

        int32_t tag() const
        {
            return detail::dense_tag;
        }

        double operator()( const point< NDim >& p ) const
        {
            index< NDim > i;
            point< NDim > s;

            // #pragma omp parallel for shared( i, s, p, cmin, hinv )
            for( unsigned d = 0; d < NDim; ++d )
            {
                double u = ( p[d] - cmin[d] ) * hinv[d];
                i[d]     = floor( u ) - 1;
                s[d]     = u - floor( u );
            }

            double f = 0;
            for( grid_iterator< NDim > d( 4 ); d; ++d )
            {
                double w = 1.0;
                for( unsigned k = 0; k < NDim; ++k )
                    w *= Bspline( d[k], s[k] );

                f += w * phi( i + ( *d ) );
            }

            return f;
        }

        void report( std::ostream& os ) const
        {
            std::ios_base::fmtflags ff( os.flags() );
            auto fp = os.precision();

            os << "dense  [" << grid[0];
            for( unsigned i = 1; i < NDim; ++i )
                os << ", " << grid[i];
            os << "] (" << phi.size() * sizeof( double ) << " bytes)";

            os.flags( ff );
            os.precision( fp );
        }

        void append_refined( const control_lattice_dense& r )
        {
            static const std::array< double, 5 > s = { 0.125, 0.500, 0.750, 0.500, 0.125 };

            for( grid_iterator< NDim > i( r.grid ); i; ++i )
            {
                double f = r.phi( *i );

                if( f == 0.0 ) continue;

                for( grid_iterator< NDim > d( 5 ); d; ++d )
                {
                    index< NDim > j;
                    bool skip = false;
                    for( unsigned k = 0; k < NDim; ++k )
                    {
                        j[k] = 2 * i[k] + d[k] - 3;
                        if( j[k] >= grid[k] )
                        {
                            skip = true;
                            break;
                        }
                    }

                    if( skip ) continue;

                    double c = 1.0;
                    for( unsigned k = 0; k < NDim; ++k )
                        c *= s[d[k]];

                    phi( j ) += f * c;
                }
            }
        }

        double fill_ratio() const
        {
            size_t total    = phi.size();
            size_t nonzeros = total - std::count( phi.data(), phi.data() + total, 0.0 );

            return static_cast< double >( nonzeros ) / total;
        }

      private:
        point< NDim > cmin, cmax, hinv;
        index< NDim > grid;

        multi_array< double, NDim > phi;
    };

    template < unsigned NDim >
    class control_lattice_sparse : public control_lattice< NDim >
    {
      public:
        template < class CooIter, class ValIter >
        control_lattice_sparse( const point< NDim >& coo_min,
                                const point< NDim >& coo_max,
                                index< NDim > grid_size,
                                CooIter coo_begin,
                                CooIter coo_end,
                                ValIter val_begin )
            : cmin( coo_min ), cmax( coo_max ), grid( grid_size )
        {
            for( unsigned i = 0; i < NDim; ++i )
            {
                hinv[i] = ( grid[i] - 1 ) / ( cmax[i] - cmin[i] );
                cmin[i] -= 1 / hinv[i];
                grid[i] += 2;
            }

            std::map< index< NDim >, two_doubles > dw;

            CooIter p = coo_begin;
            ValIter v = val_begin;

            for( ; p != coo_end; ++p, ++v )
            {
                if( !boxed( coo_min, *p, coo_max ) ) continue;

                index< NDim > i;
                point< NDim > s;

                for( unsigned d = 0; d < NDim; ++d )
                {
                    double u = ( ( *p )[d] - cmin[d] ) * hinv[d];
                    i[d]     = floor( u ) - 1;
                    s[d]     = u - floor( u );
                }

                std::array< double, power< 4, NDim >::value > w;
                double sum_w2 = 0.0;

                for( grid_iterator< NDim > d( 4 ); d; ++d )
                {
                    double prod = 1.0;
                    for( unsigned k = 0; k < NDim; ++k )
                        prod *= Bspline( d[k], s[k] );

                    w[d.position()] = prod;
                    sum_w2 += prod * prod;
                }

                for( grid_iterator< NDim > d( 4 ); d; ++d )
                {
                    double w1  = w[d.position()];
                    double w2  = w1 * w1;
                    double phi = ( *v ) * w1 / sum_w2;

                    two_doubles delta_omega = { w2 * phi, w2 };

                    append( dw[i + ( *d )], delta_omega );
                }
            }

            phi.insert( detail::make_transform_iterator( dw.begin(), delta_over_omega ),
                        detail::make_transform_iterator( dw.end(), delta_over_omega ) );
        }

        control_lattice_sparse( std::ifstream& f )
        {
            precondition( static_cast< bool >( f.read( (char*)cmin.data(), sizeof( cmin ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)cmax.data(), sizeof( cmax ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)hinv.data(), sizeof( hinv ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.read( (char*)grid.data(), sizeof( grid ) ) ), "file i/o error" );

            size_t n;
            typename sparse_grid::value_type v;

            precondition( static_cast< bool >( f.read( (char*)&n, sizeof( n ) ) ), "file i/o error" );
            for( size_t i = 0; i < n; ++i )
            {
                precondition( static_cast< bool >( f.read( (char*)&v, sizeof( v ) ) ), "file i/o error" );
                phi.insert( phi.end(), v );
            }
        }

        void write( std::ofstream& f ) const
        {
            precondition( static_cast< bool >( f.write( (char*)cmin.data(), sizeof( cmin ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)cmax.data(), sizeof( cmax ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)hinv.data(), sizeof( hinv ) ) ), "file i/o error" );
            precondition( static_cast< bool >( f.write( (char*)grid.data(), sizeof( grid ) ) ), "file i/o error" );

            size_t n = phi.size();
            precondition( static_cast< bool >( f.write( (char*)&n, sizeof( n ) ) ), "file i/o error" );
            for( auto& p : phi )
            {
                precondition( static_cast< bool >( f.write( (char*)&p, sizeof( p ) ) ), "file i/o error" );
            }
        }

        int32_t tag() const
        {
            return detail::sparse_tag;
        }

        double operator()( const point< NDim >& p ) const
        {
            index< NDim > i;
            point< NDim > s;

            for( unsigned d = 0; d < NDim; ++d )
            {
                double u = ( p[d] - cmin[d] ) * hinv[d];
                i[d]     = floor( u ) - 1;
                s[d]     = u - floor( u );
            }

            double f = 0;
            for( grid_iterator< NDim > d( 4 ); d; ++d )
            {
                double w = 1.0;
                for( unsigned k = 0; k < NDim; ++k )
                    w *= Bspline( d[k], s[k] );

                f += w * get_phi( i + ( *d ) );
            }

            return f;
        }

        void report( std::ostream& os ) const
        {
            std::ios_base::fmtflags ff( os.flags() );
            auto fp = os.precision();

            size_t grid_size = grid[0];

            os << "sparse [" << grid[0];
            for( unsigned i = 1; i < NDim; ++i )
            {
                os << ", " << grid[i];
                grid_size *= grid[i];
            }

            size_t bytes       = phi.size() * sizeof( std::pair< index< NDim >, double > );
            size_t dense_bytes = grid_size * sizeof( double );

            double compression = static_cast< double >( bytes ) / dense_bytes;
            os << "] (" << bytes << " bytes, compression: " << std::fixed << std::setprecision( 2 ) << compression
               << "), DoFs = " << phi.size();

            os.flags( ff );
            os.precision( fp );
        }

      private:
        point< NDim > cmin, cmax, hinv;
        index< NDim > grid;

        typedef std::map< index< NDim >, double > sparse_grid;
        sparse_grid phi;

        typedef std::array< double, 2 > two_doubles;

        static std::pair< index< NDim >, double > delta_over_omega( const std::pair< index< NDim >, two_doubles >& dw )
        {
            return std::make_pair( dw.first, safe_divide( dw.second[0], dw.second[1] ) );
        }

        static void append( two_doubles& a, const two_doubles& b )
        {
            a[0] += b[0];
            a[1] += b[1];
        }

        double get_phi( const index< NDim >& i ) const
        {
            typename sparse_grid::const_iterator c = phi.find( i );
            return c == phi.end() ? 0.0 : c->second;
        }
    };

}  // namespace detail

template < unsigned NDim >
class linear_approximation
{
  public:
    template < class CooIter, class ValIter >
    linear_approximation( CooIter coo_begin, CooIter coo_end, ValIter val_begin )
    {
        size_t n = std::distance( coo_begin, coo_end );
        for( auto& v : C )
            v = 0.0;

        if( n <= NDim )
        {
            // Not enough points to get a unique plane
            C[NDim] = std::accumulate( val_begin, val_begin + n, 0.0 ) / n;
            return;
        }

        double A[NDim + 1][NDim + 1] = { { 0 } };

        CooIter p = coo_begin;
        ValIter v = val_begin;

        double sum_val = 0.0;

        // Solve least-squares problem to get approximation with a plane.
        for( ; p != coo_end; ++p, ++v, ++n )
        {
            double x[NDim + 1];
            for( unsigned i = 0; i < NDim; ++i )
                x[i] = ( *p )[i];
            x[NDim] = 1.0;

            for( unsigned i = 0; i <= NDim; ++i )
            {
                for( unsigned j = 0; j <= NDim; ++j )
                {
                    A[i][j] += x[i] * x[j];
                }
                C[i] += x[i] * ( *v );
            }

            sum_val += ( *v );
        }

        // Perform LU-factorization of A in-place
        for( unsigned k = 0; k <= NDim; ++k )
        {
            if( A[k][k] == 0 ) goto singular;

            double d = 1 / A[k][k];
            for( unsigned i = k + 1; i <= NDim; ++i )
            {
                A[i][k] *= d;
                for( unsigned j = k + 1; j <= NDim; ++j )
                    A[i][j] -= A[i][k] * A[k][j];
            }
            A[k][k] = d;
        }

        for( unsigned i = 0; i <= NDim; ++i )
            if( A[i][i] == 0.0 ) goto singular;

        // Lower triangular solve:
        for( unsigned i = 0; i <= NDim; ++i )
        {
            for( unsigned j = 0; j < i; ++j )
                C[i] -= A[i][j] * C[j];
        }

        // Upper triangular solve:
        for( unsigned i = NDim + 1; i-- > 0; )
        {
            for( unsigned j = i + 1; j <= NDim; ++j )
                C[i] -= A[i][j] * C[j];
            C[i] *= A[i][i];
        }

        return;
    singular:
        for( unsigned i = 0; i < NDim; ++i )
            C[i] = 0.0;
        C[NDim] = sum_val / n;
    }

    double operator()( const point< NDim >& p ) const
    {
        double f = C[NDim];

        for( unsigned i = 0; i < NDim; ++i )
            f += C[i] * p[i];

        return f;
    }

  private:
    std::array< double, NDim + 1 > C;
};

template < unsigned NDim >
class MBA
{
  public:
    bool verbose;
    template < class CooIter, class ValIter >
    MBA( const point< NDim >& coo_min,
         const point< NDim >& coo_max,
         index< NDim > grid,
         CooIter coo_begin,
         CooIter coo_end,
         ValIter val_begin,
         unsigned max_levels                              = 8,
         double tol                                       = 1e-8,
         double min_fill                                  = 0.5,
         std::function< double( point< NDim > ) > initial = std::function< double( point< NDim > ) >() )
        : verbose( true )
    {
        init( coo_min, coo_max, grid, coo_begin, coo_end, val_begin, max_levels, tol, min_fill, initial );
    }

    template < class CooRange, class ValRange >
    MBA( const point< NDim >& coo_min,
         const point< NDim >& coo_max,
         index< NDim > grid,
         CooRange coo,
         ValRange val,
         unsigned max_levels                              = 8,
         double tol                                       = 1e-8,
         double min_fill                                  = 0.5,
         std::function< double( point< NDim > ) > initial = std::function< double( point< NDim > ) >() )
        : verbose( true )
    {
        init( coo_min, coo_max, grid, std::begin( coo ), std::end( coo ), std::begin( val ), max_levels, tol, min_fill,
              initial );
    }

    MBA( std::ifstream& f,
         std::function< double( point< NDim > ) > initial = std::function< double( point< NDim > ) >() )
        : verbose( true )
    {
        size_t n;
        detail::precondition( static_cast< bool >( f.read( (char*)&n, sizeof( n ) ) ), "file i/o error" );
        for( size_t i = 0; i < n; ++i )
        {
            int32_t tag;
            detail::precondition( static_cast< bool >( f.read( (char*)&tag, sizeof( tag ) ) ), "file i/o error" );

            switch( tag )
            {
                case detail::initial_tag:
                    detail::precondition( static_cast< bool >( initial ),
                                          "initial function definition should be provided" );
                    cl.push_back( std::make_shared< initial_approximation >( initial ) );
                    break;
                case detail::dense_tag:
                    cl.push_back( std::make_shared< dense_lattice >( f ) );
                    break;
                case detail::sparse_tag:
                    cl.push_back( std::make_shared< sparse_lattice >( f ) );
                    break;
                default:
                    detail::precondition( false, "unknown lattice tag in input file" );
            }
        }
    }

    double operator()( const point< NDim >& p ) const
    {
        double f = 0.0;

        for( auto& psi : cl )
        {
            f += ( *psi )( p );
        }

        return f;
    }

    friend std::ostream& operator<<( std::ostream& os, const MBA& h )
    {
        size_t level = 0;
        for( auto& psi : h.cl )
        {
            os << "level " << ++level << ": ";
            psi->report( os );
            os << std::endl;
        }
        return os;
    }

    void write( std::ofstream& f ) const
    {
        size_t n = cl.size();
        detail::precondition( static_cast< bool >( f.write( (char*)&n, sizeof( n ) ) ), "file i/o error" );
        for( auto& l : cl )
        {
            auto tag = l->tag();
            detail::precondition( static_cast< bool >( f.write( (char*)&tag, sizeof( tag ) ) ), "file i/o error" );
            l->write( f );
        }
    }

  private:
    typedef detail::control_lattice< NDim > lattice;
    typedef detail::initial_approximation< NDim > initial_approximation;
    typedef detail::control_lattice_dense< NDim > dense_lattice;
    typedef detail::control_lattice_sparse< NDim > sparse_lattice;

    std::list< std::shared_ptr< lattice > > cl;

    template < class CooIter, class ValIter >
    void init( const point< NDim >& cmin,
               const point< NDim >& cmax,
               index< NDim > grid,
               CooIter coo_begin,
               CooIter coo_end,
               ValIter val_begin,
               unsigned max_levels,
               double tol,
               double min_fill,
               std::function< double( point< NDim > ) > initial )
    {
        using namespace mba::detail;

        for( size_t i = 0; i < NDim; ++i )
            precondition( grid[i] > 1, "MBA: grid size in each dimension should be more than 1" );

        const ptrdiff_t n = std::distance( coo_begin, coo_end );
        std::vector< double > val( val_begin, val_begin + n );

        double res, eps = 0.0;
        for( ptrdiff_t i = 0; i < n; ++i )
            eps = std::max( eps, std::abs( val[i] ) );
        eps *= tol;

        if( initial )
        {
            // Start with the given approximation.
            cl.push_back( std::make_shared< initial_approximation >( initial ) );
            res = cl.back()->residual( coo_begin, coo_end, val.begin() );
            if( res <= eps ) return;
        }

        size_t lev = 1;
        // Create dense head of the hierarchy.
        {
            auto psi = std::make_shared< dense_lattice >( cmin, cmax, grid, coo_begin, coo_end, val.begin() );

            res         = psi->residual( coo_begin, coo_end, val.begin() );
            double fill = psi->fill_ratio();

            // print some information about the first level
            std::cout << "\nDense Level: 0\n\t";
            psi->report( std::cout );
            std::cout << ", residual = " << res;

            for( ; ( lev < max_levels ) && ( res > eps ) && ( fill > min_fill ); ++lev )
            {
                grid = grid * 2ul - 1ul;

                auto f = std::make_shared< dense_lattice >( cmin, cmax, grid, coo_begin, coo_end, val.begin() );

                res  = f->residual( coo_begin, coo_end, val.begin() );
                fill = f->fill_ratio();

                // print some information about the each level now
                std::cout << "\nDense Level: 0\n\t";
                psi->report( std::cout );
                std::cout << ", residual = " << res;

                f->append_refined( *psi );
                psi.swap( f );
            }

            cl.push_back( psi );
        }

        // Create sparse tail of the hierrchy.
        for( ; ( lev < max_levels ) && ( res > eps ); ++lev )
        {
            grid = grid * 2ul - 1ul;

            cl.push_back( std::make_shared< sparse_lattice >( cmin, cmax, grid, coo_begin, coo_end, val.begin() ) );

            res = cl.back()->residual( coo_begin, coo_end, val.begin() );

            // print information about the sparse levels
            std::cout << "\nSparse Level: " << lev << std::endl << "\t";
            cl.back()->report( std::cout );
            std::cout << ", residual = " << res;
        }
        std::cout << std::endl;
    }
};

}  // namespace mba

template < int dimension >
moab::ErrorCode compute_mba_set( moab::Interface* mbi,
                                 moab::TempestRemapper& remapper,
                                 moab::EntityHandle& meshset_source,
                                 moab::EntityHandle& meshset_target,
                                 moab::Tag& srctag,
                                 std::string target_tagname,
                                 moab::Tag& tgttag,
                                 bool normalize = true );

template <>
moab::ErrorCode compute_mba_set< 3 >( moab::Interface* mbi,
                                      moab::TempestRemapper& remapper,
                                      moab::EntityHandle& meshset_source,
                                      moab::EntityHandle& meshset_target,
                                      moab::Tag& srctag,
                                      std::string target_tagname,
                                      moab::Tag& tgttag,
                                      bool normalize )
{
    const int dimension = 3;
    moab::ErrorCode err;
    assert( target_tagname.size() );

    {
        // err = mbi->tag_get_handle( source_tagname.c_str(), 1, moab::MB_TYPE_DOUBLE, srctag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
        err = mbi->tag_get_handle( target_tagname.c_str(), 1, moab::MB_TYPE_DOUBLE, tgttag,
                                   moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( err );
    }

    moab::Range source_elems;
    err = mbi->get_entities_by_dimension( meshset_source, 2, source_elems, true );MB_CHK_ERR( err );

    moab::Range target_elems;
    err = mbi->get_entities_by_dimension( meshset_target, 2, target_elems, true );MB_CHK_ERR( err );

    const size_t nd = source_elems.size();
    const size_t ni = target_elems.size();

    std::vector< double > source_xyz, source_tdata, target_xyz, target_tdata;
    source_xyz.resize( nd * dimension, 0.0 );
    err = mbi->get_coords( source_elems, source_xyz.data() );MB_CHK_ERR( err );

    source_tdata.resize( nd, 0.0 );
    err = mbi->tag_get_data( srctag, source_elems, source_tdata.data() );MB_CHK_ERR( err );

    target_xyz.resize( ni * dimension, 0.0 );
    err = mbi->get_coords( target_elems, target_xyz.data() );MB_CHK_ERR( err );

    // Loop over all Faces in meshOverlap
    double dTotalFieldIntegralIn = 0.0, dTotalFieldIntegralOut = 0.0;
    Mesh& meshOverlap = *remapper.GetMesh( moab::Remapper::OverlapMesh );
    // Mesh& meshInput   = *remapper.GetMesh( moab::Remapper::CoveringMesh );
    // Mesh& meshOutput  = *remapper.GetMesh( moab::Remapper::TargetMesh );
    if( normalize )
    {
        // Compute intersections with MOAB with either the Kd-tree or the advancing front algorithm
        // std::cout << "Setup and compute mesh intersections between source (MPAS) and target (ROMS) meshes\n";
        // err = remapper.ComputeOverlapMesh( true, false );MB_CHK_ERR( err );

        // meshInput.ConstructEdgeMap();
        // meshOutput.ConstructEdgeMap();
        // bool concaveMeshA = false, concaveMeshB = false, allowNoOverlap = true, verbose = false;
        // int err = GenerateOverlapWithMeshes( meshInput, meshOutput, meshOverlap, "" /*outFilename*/, "Netcdf4", "exact",
        //                                      concaveMeshA, concaveMeshB, allowNoOverlap, verbose );
        // if( err )
        // {
        //     MB_CHK_SET_ERR( MB_FAILURE, "TempestRemap: Can't compute the intersection of meshes on the sphere" );
        // }

        // EntityHandle& intersection_set = remapper.GetMeshSet( moab::Remapper::OverlapMesh );

        // Loop through all overlap faces associated with this source face
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iSourceFace = meshOverlap.vecSourceFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iSourceFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralIn += source_tdata[iSourceFace] * meshOverlap.vecFaceArea[j];
        }
    }

    target_tdata.resize( ni, 0.0 );
    std::cout << "Computing the MBA interpolant now";
    // err = compute_mba( source_xyz, source_tdata, target_xyz, target_tdata );MB_CHK_ERR( err );

    // now set the data on target instance of MOAB tag
    {
        // Bounding box containing the data points.
        mba::point< dimension > lo = { -1, -1, -1 };
        mba::point< dimension > hi = { 1, 1, 1 };

        // Initial grid size.
        // const size_t init_grid_size = static_cast< size_t >( std::max( 10.0, std::sqrt( nd ) / 25 ) );
        // mba::index< 3 > grid        = { init_grid_size, init_grid_size, init_grid_size };
        mba::index< dimension > grid = { 5, 5, 2 };

        std::vector< mba::point< 3 > > coords( nd );
        size_t offset = 0;
        for( size_t k = 0; k < nd; k++, offset += 3 )
            coords[k] = mba::point< dimension >{ source_xyz[offset], source_xyz[offset + 1], source_xyz[offset + 2] };

        // Algorithm setup.
        mba::MBA< dimension > interp( lo, hi, grid, coords, source_tdata, 10 /*levels*/, 1e-14 /*tolerance*/,
                                      0.5 /*min_fill*/ );
        // mba::linear_approximation< 3 > interp( coords.begin(), coords.end(), fd.begin() );

        // Get interpolated value at arbitrary location.
        offset = 0;
        for( size_t k = 0; k < ni; k++, offset += dimension )
            target_tdata[k] =
                interp( mba::point< dimension >{ target_xyz[offset], target_xyz[offset + 1], target_xyz[offset + 2] } );
    }

    if( normalize )
    {
        // Loop through all overlap-target faces and compute integral
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iTargetFace = meshOverlap.vecTargetFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iTargetFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralOut += target_tdata[iTargetFace] * meshOverlap.vecFaceArea[j];
        }

        if( dTotalFieldIntegralOut > 1e-8 )
        {
            const double normFactor = dTotalFieldIntegralIn / dTotalFieldIntegralOut;
            for( size_t ind = 0; ind < target_tdata.size(); ind++ )
                target_tdata[ind] *= normFactor;
        }
    }

    // now set the final projected tag value on target mesh
    err = mbi->tag_set_data( tgttag, target_elems, target_tdata.data() );MB_CHK_ERR( err );

    return moab::MB_SUCCESS;
}

void spherical_to_cart( double lat, double lon, double R, double res[3] )
{
    lat *= 3.14159265358979323846 / 180;
    lon *= 3.14159265358979323846 / 180;
    res[0] = R * cos( lat ) * cos( lon );  // x coordinate
    res[1] = R * cos( lat ) * sin( lon );  // y
    res[2] = R * sin( lat );               // z
}

double mlat[2] = { 4 * M_PI, -4 * M_PI }, mlon[2] = { 4 * M_PI, -4 * M_PI };

void cart_to_spherical( const double cart3d[3], mba::point< 2 >& latlon )
{
    double c3d[3] = { cart3d[0], cart3d[1], cart3d[2] };
    moab::IntxUtils::transform_coordinates( c3d, 2 /* int projection_type */ );

    latlon[0] = c3d[1];
    latlon[1] = c3d[0];

    return;
}

template <>
moab::ErrorCode compute_mba_set< 2 >( moab::Interface* mbi,
                                      moab::TempestRemapper& remapper,
                                      moab::EntityHandle& meshset_source,
                                      moab::EntityHandle& meshset_target,
                                      moab::Tag& srctag,
                                      std::string target_tagname,
                                      moab::Tag& tgttag,
                                      bool normalize )
{
    const int dimension = 2;
    moab::ErrorCode err;
    assert( target_tagname.size() );

    {
        // err = mbi->tag_get_handle( source_tagname.c_str(), 1, moab::MB_TYPE_DOUBLE, srctag, moab::MB_TAG_DENSE );MB_CHK_ERR( err );
        err = mbi->tag_get_handle( target_tagname.c_str(), 1, moab::MB_TYPE_DOUBLE, tgttag,
                                   moab::MB_TAG_DENSE | moab::MB_TAG_CREAT );MB_CHK_ERR( err );
    }

    moab::Range source_elems;
    err = mbi->get_entities_by_dimension( meshset_source, 2, source_elems, true );MB_CHK_ERR( err );

    moab::Range target_elems;
    err = mbi->get_entities_by_dimension( meshset_target, 2, target_elems, true );MB_CHK_ERR( err );

    const size_t nd = source_elems.size();
    const size_t ni = target_elems.size();

    std::vector< double > source_xyz, source_tdata, target_xyz, target_tdata;
    source_xyz.resize( nd * 3, 0.0 );
    err = mbi->get_coords( source_elems, source_xyz.data() );MB_CHK_ERR( err );

    source_tdata.resize( nd, 0.0 );
    err = mbi->tag_get_data( srctag, source_elems, source_tdata.data() );MB_CHK_ERR( err );

    target_xyz.resize( nd * 3, 0.0 );
    err = mbi->get_coords( target_elems, target_xyz.data() );MB_CHK_ERR( err );

    // Loop over all Faces in meshOverlap
    double dTotalFieldIntegralIn = 0.0, dTotalFieldIntegralOut = 0.0;
    Mesh& meshOverlap = *remapper.GetMesh( moab::Remapper::OverlapMesh );
    if( normalize )
    {
        // Loop through all overlap faces associated with this source face
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iSourceFace = meshOverlap.vecSourceFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iSourceFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralIn += source_tdata[iSourceFace] * meshOverlap.vecFaceArea[j];
        }
    }

    target_tdata.resize( ni, 0.0 );
    std::cout << "Computing the MBA interpolant now";
    // now set the data on target instance of MOAB tag
    {
        std::vector< mba::point< dimension > > coords( nd );
        size_t offset = 0;
        for( size_t k = 0; k < nd; k++, offset += 3 )
            cart_to_spherical( &source_xyz[offset], coords[k] );

        // const double dPadding = M_PI;
        // Bounding box containing the data points.
        // mba::point< dimension > lo = { mlon[0] - dPadding, mlat[0] - dPadding };
        // mba::point< dimension > hi = { mlon[1] + dPadding, mlat[1] + dPadding };
        mba::point< dimension > lo = { -4 * M_PI, -4 * M_PI };
        mba::point< dimension > hi = { 4 * M_PI, 4 * M_PI };

        // Initial grid size.
        // const size_t init_grid_size = static_cast< size_t >( std::max( 10.0, std::sqrt( nd ) / 10 ) );
        // mba::index< 2 > grid        = { init_grid_size, init_grid_size };
        mba::index< dimension > grid = { 5, 5 };

        std::cout << "Minimum latitude = " << mlat[0] << " and longitude = " << mlon[0] << std::endl;
        std::cout << "Maximum latitude = " << mlat[1] << " and longitude = " << mlon[1] << std::endl;

        // Algorithm setup.
        mba::MBA< dimension > interp( lo, hi, grid, coords, source_tdata, 10 /*levels*/, 1e-12 /*tolerance*/,
                                      0.5 /*min_fill*/ );
        // mba::linear_approximation< 2 > interp( coords.begin(), coords.end(), fd.begin() );

        // Get interpolated value at arbitrary location.
        offset = 0;
        mba::point< dimension > icoords;
        for( size_t k = 0; k < ni; k++, offset += 3 )
        {
            cart_to_spherical( &target_xyz[offset], icoords );
            target_tdata[k] = interp( icoords );
        }
    }

    if( normalize )
    {
        // Loop through all overlap-target faces and compute integral
        for( size_t j = 0; j < meshOverlap.faces.size(); j++ )
        {
            int iTargetFace = meshOverlap.vecTargetFaceIx[j];

            // signal to not participate, because it is a ghost target
            if( iTargetFace < 0 ) continue;  // skip and do not do anything

            dTotalFieldIntegralOut += target_tdata[iTargetFace] * meshOverlap.vecFaceArea[j];
        }

        if( dTotalFieldIntegralOut > 1e-8 )
        {
            const double normFactor = dTotalFieldIntegralIn / dTotalFieldIntegralOut;
            for( size_t ind = 0; ind < target_tdata.size(); ind++ )
                target_tdata[ind] *= normFactor;
        }
    }

    // now set the final projected tag value on target mesh
    err = mbi->tag_set_data( tgttag, target_elems, target_tdata.data() );MB_CHK_ERR( err );

    return moab::MB_SUCCESS;
}

///////////////////////////////////////////////
#endif
