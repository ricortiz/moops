#ifndef KRYLOV_MATH_HPP
#define KRYLOV_MATH_HPP

/**
 * \brief Returns the values of the rotation matrix that perform a Givens' rotation.
 *   That is, it returns the following transformation
 *   1/sqrt(dx^2+dy^2)[dx dy] [dx] = [dx^2+dy^2]/sqrt(dx^2+dy^2)
 *   1/sqrt(dx^2+dy^2)[-dy dx][dy]   [0]
 *
 *          /|
 *        /  | dy
 *      /____|
 *        dx
 *
 * \param dx The length of triangle in x.  cos(theta) = dx/sqrt(dx^2+dy^2)
 * \param dy This value will transform to zero. sin(theta) = dy/sqrt(dx^2+dy^2)
 * \param cs dx/sqrt(dx^2+dy^2)
 * \param sn dy/sqrt(dx^2+dy^2)
 *
 **/
template<typename value_type>
inline void get_rotation ( const value_type &dx, const value_type &dy, value_type &cs, value_type &sn )
{
    if ( dy == 0.0 )
    {
        cs = 1.0;
        sn = 0.0;
    }
    else if ( std::fabs ( dy ) > std::fabs ( dx ) )
    {
        assert ( dy != 0. );
        value_type temp = dx / dy;
        sn = 1.0 / std::sqrt ( 1.0 + temp * temp );
        cs = temp * sn;
    }
    else
    {
        assert ( dx != 0. );
        value_type temp = dy / dx;
        cs = 1.0 / std::sqrt ( 1.0 + temp * temp );
        sn = temp * cs;
    }
}

/**
 * \brief Apply the Givens rotation.  Does the rotation matrix-vector product.
 *
 * \param dx Result of the rotation
 * \param dy Result of the rotation, should equal 0.
 * \param cs Rotation cosine parameter
 * \param sn Rotation sine parameter
 *
 **/
template<typename value_type>
inline void apply_rotation ( value_type &dx, value_type &dy, const value_type &cs, const value_type &sn )
{
    value_type temp  =  cs * dx + sn * dy;
    dy = cs * dy - sn * dx;
    dx = temp;
}

/**
 * \brief Back substitution solver for an upper triangular linear system. Uy=b.
 *  Where U = H(0:k,0:k)  and b = g(0:k).
 *
 * \param k Index in the GMRES loop
 * \return Upon return, g(0:k) contains the solution.
 *
 **/
inline void back_solve ( int k )
{
    for ( int i = k-1; i >= 0; i-- )
    {
        derived().g ( i ) /= derived().H ( i,i );
        for ( int j = i - 1; j >= 0; j-- )
            derived().g ( j ) -= derived().H ( j,i ) * derived().g ( i );
    }
}

template<typename value_type>
inline value_type dot ( const value_type *x, const value_type *y )
{
    value_type s = value_type ( 0 );
    for ( size_t i = 0; i < m_size; ++i )
        s += x[i]*y[i];
    return s;
}

inline value_type norm ( const value_type *x )
{
    return std::sqrt ( dot ( x,x ) );
}

template<typename value_type, int ACCUM_N> __device__ __host__
value_type dot ( const value_type* x, const value_type *y, size_t size )
{
    //Accumulators cache
    __shared__ value_type accumResult[ACCUM_N];
    
    ////////////////////////////////////////////////////////////////////////
    // Each accumulator cycles through vectors with
    // stride equal to number of total number of accumulators ACCUM_N
    // At this stage ACCUM_N is only preferred be a multiple of warp size
    // to meet memory coalescing alignment constraints.
    ////////////////////////////////////////////////////////////////////////
    for (int iAccum = threadIdx.x; iAccum < ACCUM_N; iAccum += blockDim.x) {
        value_type sum = 0;
        
        for (int pos = iAccum; pos < size; pos += ACCUM_N)
            sum += x[pos] * y[pos];
        
        accumResult[iAccum] = sum;
    }
    
    ////////////////////////////////////////////////////////////////////////
    // Perform tree-like reduction of accumulators' results.
    // ACCUM_N has to be power of two at this stage
    ////////////////////////////////////////////////////////////////////////
    for (int stride = ACCUM_N / 2; stride > 0; stride >>= 1) {
        __syncthreads();
        for (int iAccum = threadIdx.x; iAccum < stride; iAccum += blockDim.x)
            accumResult[iAccum] += accumResult[stride + iAccum];
    }
    
    return accumResult[0];
}

template<typename value_type, int ACCUM_N> __device__ __host__
value_type norm ( const value_type* x, size_t size )
{
    return sqrt(dot<value_type,ACCUM_N>(x,x,size));    
}

#endif
