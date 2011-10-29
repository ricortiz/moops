#ifndef SDC_STORAGE_RAW_HPP
#define SDC_STORAGE_RAW_HPP

#include<algorithm>

namespace SDC
{
    enum {EXPLICIT,SEMI_IMPLICIT,MULTIRATE,SEMI_IMPLICIT_MULTIRATE};
}


/** \internal
 * Static array
 */
template<typename T, int solution_size, int function_size>
struct sdc_arrays
{
    T* X[solution_size];
    T* F[function_size];
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
template<typename T, int sdc_nodes, int multirate_nodes, int sdc_type>
class sdc_storage;

// purely fixed-size sdc arrays
template<typename T, int sdc_nodes, int multirate_nodes, int sdc_type>
class sdc_storage
{
        sdc_arrays<T,sdc_nodes,sdc_nodes> m_data;
        size_t m_ode_size;

    public:
        inline explicit sdc_storage ( size_t ode_size ) : m_ode_size ( ode_size )
        {
            for ( int i = 0; i < sdc_nodes; ++i )
            {
                m_data.X[i] = new T[ode_size];
                m_data.F[i] = new T[ode_size];
            }
        }

        ~sdc_storage()
        {
            for ( int i = 0; i < sdc_nodes; ++i )
            {
                delete [] m_data.X[i];
                delete [] m_data.F[i];
            }
        }
        inline void swap ( sdc_storage& other ) { std::swap ( m_data,other.m_data ); }
        inline void update()
        {
            std::copy ( m_data.F[sdc_nodes-1],m_data.F[sdc_nodes-1]+m_ode_size,m_data.F[0] );
            std::copy ( m_data.X[sdc_nodes-1],m_data.X[sdc_nodes-1]+m_ode_size,m_data.X[0] );
        }

        inline void init(T *x, T *Fx)
        {
            std::copy ( x,x+m_ode_size,m_data.X[0] );
            std::copy ( Fx,Fx+m_ode_size,m_data.F[0] );
        }
        
        inline const T **F() const { return m_data.F; }
        inline T **F() { return m_data.F; }
        inline const T **X() const { return m_data.X; }
        inline T **X() { return m_data.X; }
};

template<typename T, int sdc_nodes>
class sdc_storage<T,sdc_nodes,0,SDC::EXPLICIT>
{
        sdc_arrays<T,sdc_nodes,sdc_nodes> m_data;
        size_t m_ode_size;
    public:
        inline explicit sdc_storage ( size_t ode_size ) : m_ode_size ( ode_size )
        {
            for ( int i = 0; i < sdc_nodes; ++i )
            {
                m_data.X[i] = new T[ode_size];
                m_data.F[i] = new T[ode_size];
            }
        }

        ~sdc_storage()
        {
            for ( int i = 0; i < sdc_nodes; ++i )
            {
                delete [] m_data.X[i];
                delete [] m_data.F[i];
            }
        }
        inline void swap ( sdc_storage& other ) { std::swap ( m_data,other.m_data ); }
        inline void update()
        {
            std::copy ( m_data.F[sdc_nodes-1],m_data.F[sdc_nodes-1]+m_ode_size,m_data.F[0] );
            std::copy ( m_data.X[sdc_nodes-1],m_data.X[sdc_nodes-1]+m_ode_size,m_data.X[0] );
        }

        inline void init(T *x, T *Fx)
        {
            std::copy ( x,x+m_ode_size,m_data.X[0] );
            std::copy ( Fx,Fx+m_ode_size,m_data.F[0] );
        }
        
        inline const T **F() const { return m_data.F; }
        inline T **F() { return m_data.F; }
        inline const T **X() const { return m_data.X; }
        inline T **X() { return m_data.X; }
};

template<typename T, int sdc_nodes>
class sdc_storage<T,sdc_nodes,0,SDC::SEMI_IMPLICIT>
{
        sdc_arrays<T,sdc_nodes,2*sdc_nodes> m_data;
        size_t m_ode_size;
    public:
        inline explicit sdc_storage ( size_t ode_size ) : m_ode_size ( ode_size )
        {
            for ( int i = 0; i < sdc_nodes; ++i )
            {
                m_data.X[i] = new T[ode_size];
                m_data.F[i] = new T[ode_size];
                m_data.F[i+sdc_nodes] = new T[ode_size];
            }
        }

        ~sdc_storage()
        {
            for ( int i = 0; i < sdc_nodes; ++i )
            {
                delete [] m_data.X[i];
                delete [] m_data.F[i];
                delete [] m_data.F[i+sdc_nodes];
            }
        }

        inline void swap ( sdc_storage& other ) { std::swap ( m_data,other.m_data ); }

        inline void update()
        {
            std::copy ( m_data.F[sdc_nodes-1],m_data.F[sdc_nodes-1]+m_ode_size,m_data.F[0] );
            std::copy ( m_data.F[2*sdc_nodes-1],m_data.F[2*sdc_nodes-1]+m_ode_size,m_data.F[sdc_nodes] );
            std::copy ( m_data.X[sdc_nodes-1],m_data.X[sdc_nodes-1]+m_ode_size,m_data.X[0] );
        }

        inline void init(T *x, T *Fx_i, T *Fx_e)
        {
            std::copy ( x,x+m_ode_size,m_data.X[0] );
            std::copy ( Fx_i,Fx_i+m_ode_size,m_data.F[0] );
            std::copy ( Fx_e,Fx_e+m_ode_size,m_data.F[sdc_nodes] );
        }
        
        inline const T **Fi() const { return m_data.F; }
        inline T **Fi() { return m_data.F; }
        inline const T **Fe() const { return &m_data.F[sdc_nodes]; }
        inline T **Fe() { return &m_data.F[sdc_nodes]; }
        inline const T **X() const { return m_data.X; }
        inline T **X() { return m_data.X; }

};

// template<typename T, int sdc_nodes, int multirate_nodes>
// class sdc_storage<T,sdc_nodes,multirate_nodes,SDC::MULTIRATE>
// {
//         sdc_arrays<T,sdc_nodes, ( sdc_nodes-1 ) * ( multirate_nodes-1 ) + 1> m_data;
//     public:
//         inline explicit sdc_storage() {}
//
//         inline void swap ( sdc_storage& other ) { std::swap ( m_data,other.m_data ); }
//
//         inline void update() { m_data.F[0] = m_data.F[sdc_nodes-1]; m_data.X[0] = m_data.X[sdc_nodes-1]; }
//
//         inline const T *F() const { return m_data.F; }
//
//         inline T *F() { return m_data.F; }
//
//         inline const T *X() const { return m_data.X; }
//
//         inline T *X() { return m_data.X; }
// };
//
// template<typename T, int sdc_nodes, int multirate_nodes>
// class sdc_storage<T,sdc_nodes,multirate_nodes,SDC::SEMI_IMPLICIT_MULTIRATE>
// {
//         sdc_arrays<T,sdc_nodes,2* ( ( sdc_nodes-1 ) * ( multirate_nodes-1 ) + 1 ) +sdc_nodes> m_data;
//     public:
//         inline explicit sdc_storage() {}
//
//         inline void swap ( sdc_storage& other ) { std::swap ( m_data,other.m_data ); }
//
//         inline void update() { m_data.F[0] = m_data.F[sdc_nodes-1]; m_data.X[0] = m_data.X[sdc_nodes-1]; }
//
//         inline const T *F() const { return m_data.F; }
//
//         inline T *F() { return m_data.F; }
//
//         inline const T *X() const { return m_data.X; }
//
//         inline T *X() { return m_data.X; }
// };

#endif
