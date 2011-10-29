#ifndef IO_UTILS_H
#define IO_UTILS_H

#include <sstream>

namespace IO
{
    namespace detail
    {
        template<typename vector3_type>
        void str2vec ( vector3_type & to, std::string const & from )
        {
            if ( from.empty() )
                return;

            std::istringstream ist ( from );

            char dummy;

            ist >> dummy >> to ( 0 ) >> dummy >> to ( 1 ) >> dummy >> to ( 2 ) >> dummy;
        }
        
        template<typename quaternion_type>
        void str2quat( quaternion_type & to, std::string const & from )
        {
            if ( from.empty() )
                return;

            std::istringstream ist ( from );

            char dummy;

            
            ist >> dummy >> to.s() >> dummy >> to.v() ( 0 ) >> dummy >> to.v() ( 1 ) >> dummy >> to.v() ( 2 ) >> dummy;
        }

    }
}
#endif
