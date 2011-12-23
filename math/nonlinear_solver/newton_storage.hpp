#ifndef NEWTON_STORAGE_HPP
#define NEWTON_STORAGE_HPP

/** \internal
 *
 */
template<typename T>
struct newton_arrays
{
    T *dx;
    T *f;
    T *ft;
    T *xt;
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
template<typename T>
class newton_storage;

// purely fixed-size arrays
template<typename T>
class newton_storage
{
        newton_arrays<T> m_data;
    public:
        inline explicit newton_storage (size_t size)
        {
            m_data.dx = new T[size];
            m_data.f = new T[size];
            m_data.ft = new T[size];
            m_data.xt = new T[size];
            clear(size);
        }
        void clear(size_t size)
        {
            std::fill(m_data.dx,m_data.dx+size,0.0);
            std::fill(m_data.f,m_data.f+size,0.0);
            std::fill(m_data.ft,m_data.ft+size,0.0);
            std::fill(m_data.xt,m_data.xt+size,0.0);
        }
        ~newton_storage()
        {
            delete [] m_data.dx;
            delete [] m_data.f;
            delete [] m_data.ft;
            delete [] m_data.xt;
        }
        inline void swap ( newton_storage& other ) { std::swap ( m_data,other.m_data ); }

        inline T *dx() { return m_data.dx; }
        inline T &dx(size_t i) { return * ( m_data.dx+i); }
        inline T *f ( ) { return m_data.f; }
        inline T *ft ( ) { return m_data.ft; }
        inline T *xt ( ) { return m_data.xt; }
        inline T &xt ( size_t i ) { return * ( m_data.xt+i ); }


};


#endif
