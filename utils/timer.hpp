#ifndef TIMER_HPP
#define TIMER_HPP
//
// OpenTissue Template Library
// - A generic toolbox for physics-based modeling and simulation.
// Copyright (C) 2008 Department of Computer Science, University of Copenhagen.
//
// OTTL is licensed under zlib: http://opensource.org/licenses/zlib-license.php
//

#ifdef WIN32
# define NOMINMAX
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
# undef WIN32_LEAN_AND_MEAN
# undef NOMINMAX
#else
# include<sys/time.h>
#endif

#include <map>
#include <string>


/**
 * High Resoultion Timer.
 * Based on http://www-106.ibm.com/developerworks/library/l-rt1/
 *
 * RunTime: High-performance programming techniques on Linux and Windows 2000
 * Setting up timing routines
 *   by
 * Edward G. Bradford (egb@us.ibm.com)
 * Senior Programmer, IBM
 * 01 Apr 2001
 *
 * Example usage (We recommand doubles):
 *
 *  Timer<double> timer;
 *
 *  timer.start()
 *  ...
 *  timer.stop()
 *  std::cout << "It took " << timer() << " seconds to do it" << std::endl;
 */
template<typename real_type>
class Timer
{
#ifdef WIN32
    private:
        LARGE_INTEGER m_start;   ///<
        LARGE_INTEGER m_end;     ///<
        LARGE_INTEGER m_freq;    ///<
        bool m_first;            ///<
    public:
        Timer(): m_first(true) {}
    public:
        void start()
        {
            if(m_first)
            {
                QueryPerformanceFrequency(&m_freq);
                m_first = false;
            }
            QueryPerformanceCounter(&m_start);
        }
        void stop()
        {
            QueryPerformanceCounter(&m_end);
        }
        real_type operator()()const
        {
            real_type end = static_cast<real_type>(m_end.QuadPart);
            real_type start = static_cast<real_type>(m_start.QuadPart);
            real_type freq = static_cast<real_type>(m_freq.QuadPart);
            return (end - start) / freq;
        }
#else
    private:
        std::map<std::string, struct timeval>  m_start_times;
        struct timeval m_end;     ///<

    public:
        inline void start(std::string event)
        {
            gettimeofday(&m_start_times[event], NULL);
        }
        inline real_type stop(std::string event)
        {
            gettimeofday(&m_end, NULL);
            return operator()(event);
        }
        inline real_type operator()(std::string event) 
        {
            real_type t1 =  static_cast<real_type>(m_start_times[event].tv_sec) + 1e-6 * static_cast<real_type>(m_start_times[event].tv_usec);
            real_type t2 =  static_cast<real_type>(m_end.tv_sec) + 1e-6 * static_cast<real_type>(m_end.tv_usec);
            return t2 - t1;
        }
        void reset() { m_start_times.clear(); }
        void eraseEvent(std::string event) { m_start_times.erase(event); }
#endif
};

#endif

