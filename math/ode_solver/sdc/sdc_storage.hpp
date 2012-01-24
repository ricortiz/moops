#ifndef SDC_STORAGE_RAW_HPP
#define SDC_STORAGE_RAW_HPP

#include<algorithm>

namespace SDC
{
    enum sdc_type {EXPLICIT, SEMI_IMPLICIT, MULTIRATE, SEMI_IMPLICIT_MULTIRATE};
}


/** \internal
 * Static array
 */
template<typename T, int sdc_nodes, int function_size>
struct sdc_arrays
{
    T X[sdc_nodes];
    T F[function_size];
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
template<typename T, int sdc_nodes, int multirate_nodes, SDC::sdc_type sdc>
class sdc_storage;

template<typename T, int sdc_nodes>
class sdc_storage<T, sdc_nodes, 0, SDC::EXPLICIT>
{
        sdc_arrays<T*, sdc_nodes, sdc_nodes> m_data;
    public:
        size_t ode_size;
        inline explicit sdc_storage(size_t _ode_size) : ode_size(_ode_size)
        {
            for(int i = 1; i < sdc_nodes; ++i)
            {
                m_data.X[i] = new T[ode_size];
                m_data.F[i] = new T[ode_size];
            }
        }

        ~sdc_storage()
        {
            for(int i = 1; i < sdc_nodes; ++i)
            {
                delete [] m_data.X[i];
                delete [] m_data.F[i];
            }
        }
        inline void swap(sdc_storage &other) { std::swap(m_data, other.m_data); }
        inline void update()
        {
            std::copy(m_data.F[sdc_nodes - 1], m_data.F[sdc_nodes - 1] + ode_size, m_data.F[0]);
            std::copy(m_data.X[sdc_nodes - 1], m_data.X[sdc_nodes - 1] + ode_size, m_data.X[0]);
        }

        inline void init(T *x, T *Fx)
        {
            setX0(x);
            setF0(Fx);
        }
        inline void setX0(T *x)
        {
            m_data.X[0] = x;
        }
        inline void setF0(T *Fx)
        {
            m_data.F[0] = Fx;
        }

        inline const T **F() const { return m_data.F; }
        inline T **F() { return m_data.F; }
        inline const T **X() const { return m_data.X; }
        inline T **X() { return m_data.X; }
};

template<typename T, int sdc_nodes>
class sdc_storage<T, sdc_nodes, 0, SDC::SEMI_IMPLICIT>
{
        sdc_arrays<T*, sdc_nodes, 2 * sdc_nodes> m_data;
    public:
        size_t ode_size;
        inline explicit sdc_storage(size_t _ode_size) : ode_size(_ode_size)
        {
            for(int i = 1; i < sdc_nodes; ++i)
            {
                m_data.X[i] = new T[ode_size];
                m_data.F[i] = new T[ode_size];
                m_data.F[i + sdc_nodes] = new T[ode_size];
            }
        }

        ~sdc_storage()
        {
            for(int i = 1; i < sdc_nodes; ++i)
            {
                delete [] m_data.X[i];
                delete [] m_data.F[i];
                delete [] m_data.F[i + sdc_nodes];
            }
        }

        inline void swap(sdc_storage &other) { std::swap(m_data, other.m_data); }

        inline void update()
        {
            std::copy(m_data.F[sdc_nodes - 1], m_data.F[sdc_nodes - 1] + ode_size, m_data.F[0]);
            std::copy(m_data.F[2 * sdc_nodes - 1], m_data.F[2 * sdc_nodes - 1] + ode_size, m_data.F[sdc_nodes]);
            std::copy(m_data.X[sdc_nodes - 1], m_data.X[sdc_nodes - 1] + ode_size, m_data.X[0]);
        }

        inline void init(T *x, T *Fx_i, T *Fx_e)
        {

            m_data.X[0] = x;
            m_data.F[0] = Fx_i;
            m_data.F[sdc_nodes] = Fx_e;
        }
        inline void setX0(T *x)
        {
            m_data.X[0] = x;
        }
        inline void setF0(T *Fi, T *Fe)
        {
            m_data.F[0] = Fi;
	    m_data.F[sdc_nodes] = Fe;
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
//         sdc_arrays<T*,sdc_nodes, ( sdc_nodes-1 ) * ( multirate_nodes-1 ) + 1> m_data;
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
