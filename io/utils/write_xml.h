#ifndef WRITE_XML_H
#define WRITE_XML_H
/// C++ Interface: write_xml
///
/// @brief Description: Write xml formated data
///
/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2008
/// $Id $
#include <boost/cast.hpp>  // needed for boost::numeric_cast
#include <string>
#include <iostream>
#include <fstream>
namespace IO
{
    template<typename particle_iterator, typename real_type>
    bool write_xml (
        std::string const & filename
        , particle_iterator const & particle_begin
        , particle_iterator const & particle_end
        , size_t size
        , real_type const  &time_step
        , real_type const  &time
        , bool compress = false
    )
    {
        typedef size_t           size_type;

        std::ofstream xmlfile ( filename.c_str(), std::ios::out );

        if ( !xmlfile )
        {
            std::cerr << "write_xml(): Error unable to create file: " << filename << std::endl;
            return false;
        }

        ///< GStime = global simulation time
        std::string tab = "\t";

        std::string quote = "\"";

        xmlfile.precision ( 16 );

        xmlfile << "<configuration timestep=" << quote << time_step << quote
        << " GStime=" << quote << time << quote << ">" << std::endl;


        //--- Write particles
        xmlfile << tab
        << "<particles num_particles=" << quote
        << size << quote
        << ">" << std::endl;

        size_type i = 0;

        particle_iterator p = particle_begin;

        for ( ; p != particle_end; ++p )
        {
            xmlfile << tab << tab
            << "<particle position=" << quote << p->old_position() << quote
            << " velocity="      << quote << p->velocity() << quote
            << " >" << std::endl;
            xmlfile << tab << tab << tab
            << "<force>" <<  p->force()
            << "</force>" << std::endl;
            xmlfile << tab << tab
            << "</particle>" << std::endl;
        }

        xmlfile << tab

        << "</particles>" << std::endl;

        xmlfile << "</configuration>" << std::endl;
        xmlfile.flush();
        xmlfile.close();
        return true;
    }
    
    
    template<typename particle_system, typename real_type>
    bool write_xml (
        particle_system &system,
        std::string const & filename,
        real_type timestep,
        bool compress = false
    )
    {
        typedef size_t           size_type;
        typedef typename particle_system::particle_iterator particle_iterator;
        std::ofstream xmlfile ( filename.c_str(), std::ios::out );

        if ( !xmlfile )
        {
            std::cerr << "write_xml(): Error unable to create file: " << filename << std::endl;
            return false;
        }

        ///< GStime = global simulation time
        std::string tab = "\t";
        std::string quote = "\"";

        xmlfile.precision ( 16 );

        xmlfile << "<configuration timestep=" << quote << timestep << quote
                << " GStime=" << quote << system.time()-timestep << quote << ">" << std::endl;

        //--- Write particles
        xmlfile << tab
                << "<particles num_particles=" << quote << system.particles_size() 
                << quote << ">" << std::endl;

        size_type i = 0;

        particle_iterator p = system.particle_begin();

        for ( ; p != system.particle_end(); ++p )
        {
            xmlfile << tab << tab
                << "<particle position=" << quote << p->old_position() << quote
                << " velocity="      << quote << p->velocity() << quote
                << " force="      << quote << p->force() << quote
                << ">" << std::endl;
            xmlfile << tab << tab
                << "</particle>" << std::endl;
        }

        xmlfile << tab << "</particles>" << std::endl;
        xmlfile << "</configuration>" << std::endl;
        xmlfile.flush();
        xmlfile.close();
        return true;
    }
    
    template<typename vector3_type, typename real_type>
    bool write_xml (
        std::vector<vector3_type> const &vectors,
        std::string const & filename,
        char type = 'v',
        real_type timestep = 0,
        real_type time = 0        
    )
    {
        typedef size_t           size_type;
        typedef typename std::vector<vector3_type>::iterator vector_iterator;
        std::ofstream xmlfile ( filename.c_str(), std::ios::out );

        if ( !xmlfile )
        {
            std::cerr << "write_xml(): Error unable to create file: " << filename << std::endl;
            return false;
        }

        ///< GStime = global simulation time
        std::string tab = "\t";
        std::string quote = "\"";

        xmlfile.precision ( 17 );

        xmlfile << "<configuration timestep=" << quote << timestep << quote
                << " GStime=" << quote << time << quote << ">" << std::endl;

        //--- Write particles
        xmlfile << tab
                << "<particles num_particles=" << quote << vectors.size() 
                << quote << ">" << std::endl;

        size_type i = 0;

        vector_iterator p = vectors.begin();
        switch(type)
        {
            case 'p':
            {
                for ( ; p != vectors.end(); ++p )
                {
                    xmlfile << tab << tab
                        << "<particle position=" << quote << p->old_position() << quote << ">" << std::endl;
                    xmlfile << tab << tab
                        << "</particle>" << std::endl;
                }  
                break;
            }  
            
            case 'v':
            {
                for ( ; p != vectors.end(); ++p )
                {
                    xmlfile << tab << tab
                        << "<particle velocity=" << quote << p->velocity() << quote << ">" << std::endl;
                    xmlfile << tab << tab
                        << "</particle>" << std::endl;
                }   
                break;
            }
            
            case 'f':
            {
                for ( ; p != vectors.end(); ++p )
                {
                    xmlfile << tab << tab
                        << "<particle force=" << quote << p->force() << quote << ">" << std::endl;
                    xmlfile << tab << tab
                        << "</particle>" << std::endl;
                }         
                break;
            }        
            default:
                std::cerr << "Wrong vector type!" << std::endl;
        
        }


        xmlfile << tab << "</particles>" << std::endl;
        xmlfile << "</configuration>" << std::endl;
        xmlfile.flush();
        xmlfile.close();
        return true;
    }
}
#endif

// kate: indent-mode cstyle; space-indent on; indent-width 4;
