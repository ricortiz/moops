#ifndef KRYLOV_STORAGE_HPP
#define KRYLOV_STORAGE_HPP
#include <stdlib.h>

/** \internal
 *
 */
template<typename T, size_t krylov_space>
struct krylov_arrays
{
    T H[krylov_space* ( krylov_space+1 ) /2];
    T C[krylov_space+1];
    T S[krylov_space+1];
    T G[krylov_space+1];
    T *V;
    T *R;
};

/** \internal
 *
 * \class sdc_storage
 *
 * \brief Stores the data of the sdc method
 *
 * This class stores the data of fixed-size sdc vectors
 *
 */
template<typename T, size_t krylov_space>
class krylov_storage;

template<typename T, size_t krylov_space>
class krylov_storage
{
        size_t m_system_size;
        krylov_arrays<T,krylov_space> m_data;
        
    public:
        inline explicit krylov_storage ( size_t system_size ) : m_system_size(system_size)
        {
            m_data.V = new T[system_size* ( krylov_space+1 ) ];
            m_data.R = new T[system_size];
            std::fill(m_data.H,m_data.H+krylov_space* ( krylov_space+1 ) /2,0.0);
            std::fill(m_data.C,m_data.C+krylov_space+1,0.0);
            std::fill(m_data.S,m_data.S+krylov_space+1,0.0);
            std::fill(m_data.G,m_data.G+krylov_space+1,0.0);
            std::fill(m_data.V,m_data.V+system_size* ( krylov_space+1 ),0.0);
            std::fill(m_data.R,m_data.R+system_size,0.0);
        }
        ~krylov_storage()
        {
            delete [] m_data.V;
            delete [] m_data.R;
        }
        inline void swap ( krylov_storage& other ) { std::swap ( m_data,other.m_data ); }
        inline T &H ( size_t i, size_t j ) { return * ( m_data.H + i*krylov_space+j-i*(i+1)/2 ); }
        inline T *H ( ) { return m_data.H; }
        inline T *residual () { return m_data.R; }
        inline T &residual ( size_t i ) { return * ( m_data.R + i ); }
        inline T *v ( size_t i ) { return m_data.V + i*m_system_size; }
        inline T &c ( size_t i ) { return * ( m_data.C + i ); }
        inline T &s ( size_t i ) { return * ( m_data.S + i ); }
        inline T &g ( size_t i ) { return * ( m_data.G + i ); }
        inline T *g (  ) { return m_data.G; }

};


#endif
