#ifndef WRITE_FILE
#define WRITE_FILE
/// C++ Interfaces for writing xml
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
    template<typename particle_iterator>
    bool write_force (
        std::string const & filename
        , particle_iterator const & particle_begin
        , particle_iterator const & particle_end
    )
    {
	std::cout << "Writting " << filename << std::endl;
        typedef size_t           size_type;

        std::ofstream file ( filename.c_str(), std::ios::out );

        if ( !file )
        {
            std::cerr << "write_file(): Error unable to create file: " << filename << std::endl;
            return false;
        }

		//--- Write particles
        particle_iterator p = particle_begin;

        for ( ; p != particle_end; ++p )
        {
            file << p->force()(0) << " "  << p->force()(1)<< " "  << p->force()(2) <<  std::endl;
        }

        file.flush();
        file.close();
	std::cout << " ...done" << std::endl;
        return true;
    }
}
#endif

// kate: indent-mode cstyle; space-indent on; indent-width 4;
