#ifndef LOGGER_HPP
#define LOGGER_HPP

#include<iostream>
#include<fstream>
#include "timer.hpp"

class Logger
{
    protected:
        typedef std::map<std::string, double>::iterator timer_iterator;
        
    private:
        Timer<double> m_timer;
        std::map<std::string, double> m_timer_log;
        std::map<std::string, double> m_error_log;

    public:
        // Start timer for event
        inline void startTimer(std::string event)
        {
            m_timer.start(event);
        }

        // Stop timer for event
        inline double stopTimer(std::string event, bool print = false)
        {
            m_timer_log[event] += m_timer.stop(event);               // Accumulate event time to timer
            if(print) std::cout << event << " : " << m_timer_log[event] << std::endl;// Print event and timer to screen
            return m_timer(event);                        // Return the event time
        }

        // Erase event in timer
        inline void eraseEvent(std::string event)
        {
            m_timer_log.erase(event);
            m_timer.eraseEvent(event);
        }

        // Clear timer
        inline void resetTimer()
        {
            m_timer.reset();
            m_timer_log.clear();
        }

        // Print timer event
        inline void printTimerEvent(std::string event)
        {
            for(timer_iterator t = m_timer_log.begin(); t != m_timer_log.end(); ++t)
                std::cout << event << " " << m_timer_log[event] << std::endl;
        }
        
        // Print timer
        inline void printTimer()
        {
            for(timer_iterator t = m_timer_log.begin(); t != m_timer_log.end(); ++t)
                std::cout << t->first << " " << t->second << std::endl;
        }

        // Write timer
        inline void writeTimer(std::string file = "time.log")
        {
            std::ofstream fout(file.c_str());
            for(timer_iterator t = m_timer_log.begin(); t != m_timer_log.end(); ++t)
                fout << t->first << " " << t->second << std::endl;
            fout.close();
        }
};


#endif
