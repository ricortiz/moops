#ifndef READ_XML_H
#define READ_XML_H
/// C++ Interfaces for reading xml
/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2008
/// $Id: read_xml.h 193 2010-03-25 23:01:55Z slukens $

#define TIXML_USE_STL

#include "utils/error.h"
#include "tinyxml.h"
#include <boost/cast.hpp> // needed for boost::numeric_cast

#include <string>

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include "utils.h"

namespace IO
{
    namespace detail
    {

        template<typename particle_type>
        bool handle_particle_tag ( particle_type * p, TiXmlHandle const & particle_tag )
        {
            if ( particle_tag.Element()->Attribute ( "position" ) )
                str2vec ( p->old_position() ,particle_tag.Element()->Attribute ( "position" ) );
//             std::cout << p->position() << std::endl;

            if ( particle_tag.Element()->Attribute ( "velocity" ) )
                str2vec ( p->velocity() ,particle_tag.Element()->Attribute ( "velocity" ) );
// std::cout << p->velocity() << std::endl;
            if ( particle_tag.Element()->Attribute ( "force" ) )
                str2vec ( p->force() ,particle_tag.Element()->Attribute ( "force" ) );
            
            
            TiXmlElement* child = particle_tag.FirstChild ( "force" ).ToElement();
            if ( child )
                str2vec ( p->force() ,child->GetText() );
            
            if ( child=child->NextSiblingElement() )
                str2vec ( p->old_position() ,child->GetText() );
// std::cout << p->force() << std::endl;
            return true;
        }



        template<typename particle_list>
        bool handle_particles_tag ( particle_list & particles, TiXmlHandle const & particles_tag )
        {
            typedef typename particle_list::particle_type particle_type;
            typedef typename particle_list::particle_iterator particle_iterator;
            double particles_size;

            if ( particles_tag.Element()->Attribute ( "num_particles" ) )
            {
                std::istringstream timestep_stream ( particles_tag.Element()->Attribute ( "num_particles" ) );
                timestep_stream >> particles_size;
            }

            TiXmlElement * particle_tag = particles_tag.FirstChild ( "particle" ).Element();
            if ( particles.particles_size() == 0 )
            {
                for ( ;particle_tag; particle_tag = particle_tag->NextSiblingElement ( "particle" ) )
                {
                    particle_iterator p = particles.create_particle ( particle_type() );
                    if ( !handle_particle_tag ( &*p, particle_tag ) )
                        return false;
                }
            }
            else
            {
                for ( particle_iterator p = particles.particle_begin();
                        particle_tag, p != particles.particle_end();
                        particle_tag = particle_tag->NextSiblingElement ( "particle" ), ++p )
                {
                    if ( !handle_particle_tag ( &*p, particle_tag ) )
                        return false;
                }
            }
            return true;
        }



        template<typename particle_list, typename real_type>
        bool handle_config_tag ( particle_list & particles,
                                 TiXmlHandle const & config_tag,
                                 real_type & time_step, real_type & time )
        {

            //--- Attributes of <configuration> tag
            if ( config_tag.Element()->Attribute ( "timestep" ) )
            {
                std::istringstream timestep_stream ( config_tag.Element()->Attribute ( "timestep" ) );
                timestep_stream >> time_step;
            }

            if ( config_tag.Element()->Attribute ( "GStime" ) )
            {
                std::istringstream GStime_stream ( config_tag.Element()->Attribute ( "GStime" ) );
                GStime_stream >> time;
            }

            if ( !config_tag.FirstChild ( "particles" ).Element() )
                throw error ( "handle_config_tag(): Only one <particles> tag is allowed" );

            if ( !handle_particles_tag ( particles, config_tag.FirstChild ( "particles" ) ) )
                return false;

            return true;
        }



        template<typename particle_list, typename real_type>
        bool handle_document (
            particle_list & particles
            , TiXmlHandle const & document
            , real_type & time_step
            , real_type & time
        )
        {

            TiXmlElement *root = document.FirstChildElement().Element() ;

            if ( root->NextSiblingElement() )
                throw error ( "handle_document(): Only one <configuration> tag is allowed" );

            return handle_config_tag ( particles, TiXmlHandle ( root ), time_step, time );

        }

    }
    template<typename particle_list, typename real_type>
    bool read_xml (
        std::string const & filename
        , particle_list &particles
        , real_type & time_step
        , real_type & time
    )
    {

#ifdef TIXML_USE_STL
        TiXmlDocument xml_document ( filename );
#else
        TiXmlDocument xml_document ( filename.c_str() );
#endif

        if ( !xml_document.LoadFile() )
            throw error ( "read_xml(...): Error " + filename + " not found!" );

        TiXmlHandle document_handle ( &xml_document );

        if ( !detail::handle_document ( particles, document_handle, time_step, time ) )
            return false;
        return true;

    }
}
#endif
// kate: indent-mode cstyle; space-indent on; indent-width 4;
