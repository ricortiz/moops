#ifndef KRYLOV_STORAGE_CUDA_HPP
#define KRYLOV_STORAGE_CUDA_HPP

#include <cuda.h>

/** \internal
 *
 * \class sdc_storage
 *
 * \brief Stores the data of the sdc method
 *
 * This class stores the data of fixed-size sdc vectors
 *
 */
template<typename T, size_t size, size_t krylov_space>
class krylov_storage;

// purely fixed-size arrays
template<typename T, size_t size, size_t krylov_space>
class krylov_storage
{
        T *H;
        T *C;
        T *S;
        T *G;
        T *V;
        T *R;
    public:
        inline explicit krylov_storage ( )
        {
            cudaMalloc((void**)&H, sizeof(T) * krylov_space* ( krylov_space+1 ) /2);
            cudaMalloc((void**)&C, sizeof(T) * ( krylov_space+1 ) );
            cudaMalloc((void**)&S, sizeof(T) * ( krylov_space+1 ) );
            cudaMalloc((void**)&G, sizeof(T) * ( krylov_space+1 ) );
            cudaMalloc((void**)&V, sizeof(T) * size* ( krylov_space+1 ) );
            cudaMalloc((void**)&R, sizeof(T) * size );
        }
        ~krylov_storage()
        {
            cudaFree((void**)H);
            cudaFree((void**)C);
            cudaFree((void**)S);
            cudaFree((void**)G);
            cudaFree((void**)V);
            cudaFree((void**)R);
        }
        __device__ __host__
        inline T &H ( size_t i, size_t j ) { return * ( H + i*krylov_space+j-i* ( i+1 ) /2 ); }
        __device__ __host__
        inline T *H ( ) { return H; }
        __device__ __host__
        inline T *residual () { return R; }
        __device__ __host__
        inline T &residual ( size_t i ) { return * ( R + i ); }
        __device__ __host__
        inline T *v ( size_t i ) { return V + i*size; }
        __device__ __host__
        inline T &c ( size_t i ) { return * ( C + i ); }
        __device__ __host__
        inline T &s ( size_t i ) { return * ( S + i ); }
        __device__ __host__
        inline T &g ( size_t i ) { return * ( G + i ); }
        __device__ __host__
        inline T *g ( ) { return G; }

};


#endif
