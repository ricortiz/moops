#ifndef READ_CONFIG_H
#define READ_CONFIG_H
/// C++ Interfaces for reading xml configuration file
/// Author: Ricardo Ortiz <ricardo.ortiz@tulane.edu>, (C) 2009
/// $Id: read_config.h 31 2009-10-14 15:56:30Z rortiz $

#define TIXML_USE_STL
#include <queue>
#include "utils.h"
#include "utils/error.h"
#include "utils/warm_start.h"
#include "tinyxml.h"
#include <boost/lexical_cast.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <OpenTissue/core/containers/mesh/common/util/mesh_add2mesh.h>
#include <OpenTissue/core/containers/mesh/common/util/mesh_deformation_modifiers.h>

#include "geometry/body_shape.h"
#include "geometry/cilia_shape.h"
#include "geometry/circular_helix_shape.h"
#include "geometry/geometry_base.h"
#include "geometry/helix_shape.h"
#include "geometry/particle_cluster.h"
#include "geometry/sine_shape.h"
#include "geometry/sphere_shape.h"
#include "geometry/torus_shape.h"

#include "motors/spring_motor.h"

namespace fs = boost::filesystem;

namespace IO
{

    namespace detail
    {

        template<typename Geometry>
        bool handle_base_geometry_config(TiXmlHandle const & geometry_tag, Geometry &geometry, std::stringstream *stream = 0)
        {
            if (geometry_tag.FirstChild().Element()->Attribute("shape"))
            {
                int type = boost::lexical_cast<int> (geometry_tag.FirstChild().Element()->Attribute("shape"));

                switch (type)
                {
                    case 0:
                    {
                        geometry.shape() = Geometry::base_type::open;

                        if (stream) *stream << "\tGeometry shape type = open" << std::endl;

                        break;
                    }

                    case 1:
                    {
                        geometry.shape() = Geometry::base_type::closed;

                        if (stream) *stream << "\tGeometry shape type = closed" << std::endl;

                        break;
                    }

                    case 2:
                    {
                        geometry.shape() = Geometry::base_type::half;

                        if (stream) *stream << "\tGeometry shape type = half" << std::endl;

                        break;
                    }

                    case 3:
                    {
                        geometry.shape() = Geometry::base_type::cluster;

                        if (stream) *stream << "\tGeometry shape type = cluster" << std::endl;

                        break;
                    }

                    case 4:
                    {
                        geometry.shape() = Geometry::base_type::line;

                        if (stream) *stream << "\tGeometry shape type = line" << std::endl;

                        break;
                    }

                    default:
                        geometry.shape() = Geometry::base_type::open;

                        if (stream) *stream << "\tGeometry shape type = open" << std::endl;
                }
            }

            if (geometry_tag.FirstChild().Element()->Attribute("connection"))
            {
                int type = boost::lexical_cast<int> (geometry_tag.FirstChild().Element()->Attribute("connection"));

                switch (type)
                {
                    case 0:
                    {
                        geometry.connection() = Geometry::base_type::springs;

                        if (stream) *stream << "\tGeometry connection type = spring" << std::endl;

                        break;
                    }

                    case 1:
                    {
                        geometry.connection() = Geometry::base_type::variable_springs;

                        if (stream) *stream << "\tGeometry connection type = variable spring (Ghost geometry driven)" << std::endl;

                        break;
                    }

                    case 2:
                    {
                        geometry.connection() = Geometry::base_type::none;

                        if (stream) *stream << "\tGeometry connection type = no connection" << std::endl;

                        break;
                    }

                    default:
                    {
                        std::cout << "Warning! handle_base_geometry_config(): unknown connection attribute." << std::endl;
                        geometry.connection() = Geometry::base_type::none;

                        if (stream) *stream << "\tGeometry connection type = no connection" << std::endl;
                    }
                }
            }
            else
            {
                geometry.connection() = Geometry::base_type::none;

                if (stream) *stream << "\tGeometry connection type = no connection" << std::endl;
            }


            if (geometry_tag.FirstChild().Element()->Attribute("solve_forces"))
            {
                geometry.solve_forces() = boost::lexical_cast<bool> (geometry_tag.FirstChild().Element()->Attribute("solve_forces"));

                if (stream) *stream << "\tSolve for forces in geometry = true" << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("position"))
            {
                str2vec(geometry.position(),geometry_tag.FirstChild().Element()->Attribute("position"));

                if (stream) *stream << "\tPosition = " << geometry.position() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("cross_sections"))
            {
                geometry.cross_sections() = boost::lexical_cast<std::size_t> (geometry_tag.FirstChild().Element()->Attribute("cross_sections"));

                if (stream) *stream << "\tCross sections = " << geometry.cross_sections() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("points_cross_sections"))
            {
                geometry.points_cross_sections() = boost::lexical_cast<std::size_t> (geometry_tag.FirstChild().Element()->Attribute("points_cross_sections"));

                if (stream) *stream << "\tPoints per cross section = " << geometry.points_cross_sections() << std::endl;
            }

            return true;
        }

        template<typename Geometry>
        bool handle_sine_geometry_config(TiXmlHandle const & geometry_tag, Geometry &geometry, std::stringstream *stream = 0)
        {
            if (stream) *stream << "Sinusoidal Geometry" << std::endl;

            handle_base_geometry_config(geometry_tag,geometry,stream);

            if (geometry_tag.FirstChild().Element()->Attribute("speed"))
            {
                geometry.speed() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("speed"));

                if (stream) *stream << "\tWave speed = " << geometry.speed() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("amplitude"))
            {
                geometry.amplitude() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("amplitude"));

                if (stream) *stream << "\tAmplitude = " << geometry.amplitude() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("pitch"))
            {
                geometry.pitch() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("pitch"));

                if (stream) *stream << "\tPitch = " << geometry.pitch() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("length"))
            {
                geometry.length() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("length"));

                if (stream) *stream << "\tLength = " << geometry.length() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("radius"))
            {
                geometry.radius() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("radius"));

                if (stream) *stream << "\tTube Radius = " << geometry.radius() << std::endl;
            }

            return true;
        }

        template<typename Geometry>
        bool handle_circular_helix_geometry_config(TiXmlHandle const & geometry_tag, Geometry &geometry, std::stringstream *stream = 0)
        {
            if (stream) *stream << "Circular Helix Geometry" << std::endl;

            handle_base_geometry_config(geometry_tag,geometry,stream);

            if (geometry_tag.FirstChild().Element()->Attribute("circle_radius"))
            {
                geometry.circle_rad() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("circle_radius"));

                if (stream) *stream << "\tCircular Radius = " << geometry.radius() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("radius"))
            {
                geometry.radius() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("radius"));

                if (stream) *stream << "\tTube Radius = " << geometry.radius() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("speed"))
            {
                geometry.speed() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("speed"));

                if (stream) *stream << "\tWave speed = " << geometry.speed() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("wavelength"))
            {
                geometry.wavelength() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("wavelength"));

                if (stream) *stream << "\tWavelength = " << geometry.wavelength() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("amplitude"))
            {
                geometry.amplitude() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("amplitude"));

                if (stream) *stream << "\tAmplitude = " << geometry.amplitude() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("pitch"))
            {
                geometry.pitch() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("pitch"));

                if (stream) *stream << "\tPitch = " << geometry.pitch() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("length"))
            {
                geometry.length() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("length"));

                if (stream) *stream << "\tLength = " << geometry.length() << std::endl;
            }

            return true;
        }

        template<typename Geometry>
        bool handle_body_geometry_config(TiXmlHandle const & geometry_tag, Geometry &geometry, std::stringstream *stream = 0)
        {
            if (stream) *stream << "Cilindrical Body Geometry" << std::endl;

            handle_base_geometry_config(geometry_tag,geometry,stream);

            if (geometry_tag.FirstChild().Element()->Attribute("length"))
            {
                geometry.length() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("length"));

                if (stream) *stream << "\tLength = " << geometry.length() << std::endl;
            }

            if (geometry_tag.FirstChild().Element()->Attribute("radius"))
            {
                geometry.radius() = boost::lexical_cast<double> (geometry_tag.FirstChild().Element()->Attribute("radius"));

                if (stream) *stream << "\tCilinder Radius = " << geometry.radius() << std::endl;
            }

            return true;
        }

        template<typename Data>
        bool handle_geometry_tag(TiXmlHandle const & geometries_tag, Data &data)
        {
            TiXmlElement * geometry_tag = geometries_tag.FirstChild("geometry").Element();

            for (;geometry_tag; geometry_tag = geometry_tag->NextSiblingElement("geometry"))
            {
                if (!geometry_tag->Attribute("type"))
                    throw error("handle_geometry_tag(): Missing geometry type attribute");

                std::string geometry_type = geometry_tag->Attribute("type");

                if (geometry_type == "body_shape")
                {
                    handle_body_geometry_config(geometry_tag,data.body_shape);
                }
                else if (geometry_type == "sine_shape")
                    handle_sine_geometry_config(geometry_tag,data.sine_shape);

//                 else if ( geometry_type == "helix_shape" )
//                     handle_geometry_config_tag ( geometry_tag,data.helix_shape );
            }

            return true;
        }

        template<typename Motor>
        void set_motor_props(TiXmlHandle const & stokes_tag, Motor &motor, std::stringstream *stream = 0)
        {
            std::string motor_type = stokes_tag.Element()->Attribute("motor");

            if (motor_type == "sheared_spring_motor")
            {
                motor.motor_type() = Motor::sheared;

                if (stream) *stream << "\t\tType = sheared" << std::endl;

                if (stokes_tag.Element()->Attribute("shear_rate"))
                {
                    motor.shear_rate() = boost::lexical_cast<double> (stokes_tag.Element()->Attribute("shear_rate"));
                }

                if (stream) *stream << "\t\tShear rate = " << motor.shear_rate() << std::endl;
            }

//                 if(motor == "data_motor") { data.surface.motor(); }
//                 else if(motor == "rotlet_motor"){}
//                 else if(motor == "sheared_spring_motor"){}
//                 else if(motor == "spring_motor"){}
//                 else if(motor == "spring_motor_from_file"){}

        }

        template<typename Data>
        bool handle_stokes_tag(TiXmlHandle const & stokes_tag, Data &data, std::stringstream *stream = 0)
        {
            if (stream) *stream << "Stokes Solver" << std::endl;

            if (stokes_tag.Element()->Attribute("motor"))
            {
                if (stream) *stream << "\tMotor" << std::endl;

                set_motor_props(stokes_tag, data.surface.motor(),stream);
            }
            else
            {
                if (stream) *stream << "\tMotor" << std::endl;

                if (stream) *stream << "\t\tWarning! Default motor attached to Stokes solver." << std::endl;
            }

            if (stokes_tag.Element()->Attribute("with_images"))
            {
                data.surface.with_images() = boost::lexical_cast<bool> (stokes_tag.Element()->Attribute("with_images"));

                if (stream) *stream << "\tWith images = " << data.surface.with_images() << std::endl;
            }
            else
            {
                if (stream) *stream << "\tWarning! Missing images parameter in configuration file." << std::endl;
            }


            if (stokes_tag.Element()->Attribute("delta"))
            {
                double delta = boost::lexical_cast<double> (stokes_tag.Element()->Attribute("delta"));
                data.surface.delta_sqr() = delta*delta;

                if (stream) *stream << "\tDelta squared = " << data.surface.delta_sqr() << std::endl;
            }
            else
            {
                if (stream) *stream << "\tWarning! Missing delta parameter in configuration file." << std::endl;
            }

            if (stokes_tag.Element()->Attribute("init_linear_solver"))
            {
                bool solve_forces = boost::lexical_cast<bool> (stokes_tag.Element()->Attribute("init_linear_solver"));

                if (solve_forces)
                {
                    data.surface.init_linear_solver();

                    if (stream) *stream << "\tLinear Solver Initialized" << std::endl;
                }
            }
            else
            {
                if (stream) *stream << "\tWarning! Missing linear solver, default attached (need to be initialized)." << std::endl;
            }

            return true;
        }

        template<typename Data>
        bool handle_config_function_tag(TiXmlHandle const & config_function_tag, Data &data)
        {
            if (config_function_tag.Element()->Attribute("id"))
            {
                data.config_function = *config_function_tag.Element()->Attribute("id");
                return true;
            }

            return false;
        }


        template<typename Data>
        bool handle_warm_starting_tag(TiXmlHandle const & warm_starting_tag, Data &data)
        {

            if (warm_starting_tag.Element()->Attribute("path"))
            {
                data.warm_starting = 1;
                data.warm_starting_file = warm_starting_tag.Element()->Attribute("path");
                return true;
            }

            return false;
        }

        template<typename Data>
        bool handle_recording_tag(TiXmlHandle const & recording_tag, Data &data, std::stringstream *stream)
        {

            if (recording_tag.Element()->Attribute("path"))
            {
                data.recording_on = 1;
                data.data_path = recording_tag.Element()->Attribute("path");

                if (stream) *stream << "Data directory = " << data.data_path << std::endl;

                return true;
            }

            return false;
        }

        template<typename Data>
        bool handle_config_tag(TiXmlHandle const & config_tag, Data &data)
        {
            //--- Attributes of <configuration> tag
            if (config_tag.Element()->Attribute("timestep"))
                data.timestep = boost::lexical_cast<double> (config_tag.Element()->Attribute("timestep"));

            if (config_tag.Element()->Attribute("save_freq"))
                data.save_freq = boost::lexical_cast<std::size_t> (config_tag.Element()->Attribute("save_freq"));

            if (!config_tag.FirstChild("warmstarting").Element())
                data.warm_starting = 0;
            else
                if (!handle_warm_starting_tag(config_tag.FirstChild("warmstarting"), data))
                    return false;

            if (!config_tag.FirstChild("recording").Element())
                data.recording_on = 0;
            else
                if (!handle_recording_tag(config_tag.FirstChild("recording"), data))
                    return false;

            if (!config_tag.FirstChild("config_function").Element())
                data.config_function = '0';
            else
                if (!handle_config_function_tag(config_tag.FirstChild("config_function"), data))
                    return false;

            if (!config_tag.FirstChild("stokes").Element())
            {
                std::cout << "Warning! Stokes parameters missing in configuration file." << std::endl;
                return true;
            }
            else
                if (!handle_stokes_tag(config_tag.FirstChild("stokes"), data))
                    return false;

            if (!config_tag.FirstChild("geometries").Element())
            {
                std::cout << "Warning! Geometries parameters missing in configuration file." << std::endl;
                return true;
            }
            else
                if (!handle_geometry_tag(config_tag.FirstChild("geometries"), data))
                    return false;

            return true;

        }

        template<typename Data>
        bool handle_document(TiXmlHandle const & document, Data &data)
        {

            TiXmlElement *root = document.FirstChildElement().Element() ;

            if (root->NextSiblingElement())
                throw error("handle_document(): Only one <configuration> tag is allowed");

            return handle_config_tag(TiXmlHandle(root), data);

        }

    }

    template<typename Data>
    bool read_config(std::string const & filename, Data &data)
    {

#ifdef TIXML_USE_STL
        TiXmlDocument xml_document(filename);
#else
        TiXmlDocument xml_document(filename.c_str());
#endif

        if (!xml_document.LoadFile())
            throw error("read_xml(...): Error " + filename + " not found!");

        TiXmlHandle document_handle(&xml_document);

        if (!detail::handle_document(document_handle, data))
            return false;

        return true;
    }


    template<typename Types, typename Data>

    class ReadConfig : public Types
    {

        protected:
            typedef typename Types::math_types  math_types;
            typedef typename Types::mesh_type   mesh_type;
            typedef GeometryBase<Types>         geometry_type;
            typedef BodyShape<Types>            body_shape_type;
            typedef CiliaShape<Types>           cilia_shape_type;
            typedef CircularHelixShape<Types>   circular_helix_type;
            typedef HelixShape<Types>           helix_shape_type;
            typedef ParticleCluster<Types>      particle_cluster_type;
            typedef SineShape<Types>            sine_shape_type;

        protected:
            typedef typename math_types::real_type              real_type;
            typedef typename math_types::vector3_type           vector3_type;
            typedef typename math_types::matrix3x3_type         matrix3x3_type;
            typedef typename math_types::quaternion_type        quaternion_type;
            typedef typename math_types::value_traits           value_traits;
            typedef typename boost::ptr_vector<geometry_type>   geometry_container;

        private:
            Data &m_data;
            TiXmlDocument m_xml_document;
            quaternion_type m_rotation;

        public:
            ReadConfig(Data &data) : m_data(data) {}

            bool operator()(const std::string &filename, std::stringstream *stream = 0)
            {
                typedef typename geometry_container::iterator iterator;

                if (!m_xml_document.LoadFile(filename))
                    throw error("read_xml(...): Error " + m_xml_document.ValueTStr() + " not found!");

                TiXmlHandle document_handle(&m_xml_document);

                TiXmlElement *root = document_handle.FirstChildElement().Element() ;

                if (root->NextSiblingElement())
                    throw error("handle_document(): Only one <configuration> tag is allowed");

                handle_config_tag(TiXmlHandle(root),stream);

                real_type pi_half = value_traits::pi_half();

                real_type two_pi = 4.0 * pi_half;

                /// Set helix shape
                mesh_type     mesh;

                size_t range_start = 0, range_end = 0;

                iterator end = m_data.geometries.end();

                for (iterator g = m_data.geometries.begin(); g != end; ++g)
                {
                    g->connect(m_data.surface);
                    g->create(mesh);
                    OpenTissue::mesh::add2mesh(mesh, m_data.mesh);
                    range_end += mesh.size_vertices();
//                     if(stream) *stream << "Geometry Range:      =\t["     << range_start << "," << range_end << "]" << std::endl;
                    g->set_range(range_start,range_end);

                    range_start = range_end;
                    mesh.clear();
                }

                /// Apply rotation if any
                matrix3x3_type rot(m_rotation);

                OpenTissue::mesh::rotate(m_data.mesh,rot);

                m_data.surface.init(m_data.mesh, true);

                /// Check if we need to do warm starting
                this->warm_start(stream);

                // Handle stokes tag after setting geometries in case we need to build matrix
                if (!TiXmlHandle(root).FirstChild("stokes").Element())
                {
                    if (stream) *stream << "Warning! Stokes parameters missing in configuration file." << std::endl;

                    return false;
                }
                else
                    if (!detail::handle_stokes_tag(TiXmlHandle(root).FirstChild("stokes"), m_data, stream))
                        return false;

//                 create_hook(m_data.geometries[0], m_data.geometries[1], m_data);
            }

            bool handle_config_tag(TiXmlHandle const & config_tag, std::stringstream *stream = 0)
            {
                //--- Attributes of <configuration> tag
                if (config_tag.Element()->Attribute("timestep"))
                {
                    m_data.timestep = boost::lexical_cast<double> (config_tag.Element()->Attribute("timestep"));

                    if (stream) *stream << "Timestep = " << m_data.timestep << std::endl;
                }

                if (config_tag.Element()->Attribute("save_freq"))
                {
                    m_data.save_freq = boost::lexical_cast<std::size_t> (config_tag.Element()->Attribute("save_freq"));

                    if (stream) *stream << "Save data frequency = " << m_data.save_freq << std::endl;
                }
                else
                {
                    m_data.save_freq = 1;

                    if (stream) *stream << "Save data frequency set to default = " << m_data.save_freq << std::endl;
                }

                if (!config_tag.FirstChild("warmstarting").Element())
                    m_data.warm_starting = 0;
                else
                    if (!detail::handle_warm_starting_tag(config_tag.FirstChild("warmstarting"), m_data))
                        return false;

                if (!config_tag.FirstChild("recording").Element())
                    m_data.recording_on = 0;
                else
                    if (!detail::handle_recording_tag(config_tag.FirstChild("recording"), m_data, stream))
                        return false;

                if (!config_tag.FirstChild("geometries").Element())
                {
                    if (stream) *stream << "Warning! Geometries parameters missing in configuration file." << std::endl;

                    return true;
                }
                else
                    if (!handle_geometry_tag(config_tag.FirstChild("geometries"),stream))
                        return false;



                return true;

            }

            bool handle_geometry_tag(TiXmlHandle const & geometries_tag,std::stringstream *stream = 0)
            {
                if (geometries_tag.Element()->Attribute("rigidity"))
                    m_data.surface.rigidty() = 2;
                else
                    m_data.surface.rigidty() = 1;

                if (geometries_tag.Element()->Attribute("rotate"))
                {
                    IO::detail::str2quat(m_rotation,geometries_tag.Element()->Attribute("rotate"));
                    real_type angle = m_rotation.s()/2;
                    m_rotation.s() = std::cos(angle);
                    m_rotation.v() *= std::sin(angle);
                }

                TiXmlElement * geometry_tag = geometries_tag.FirstChild("geometry").Element();

                for (;geometry_tag; geometry_tag = geometry_tag->NextSiblingElement("geometry"))
                {
                    if (!geometry_tag->Attribute("type"))
                        throw error("handle_geometry_tag(): Missing geometry type attribute");

                    std::string geometry_type = geometry_tag->Attribute("type");

                    if (geometry_type == "body_shape")
                    {
                        m_data.geometries.push_back(new body_shape_type());
                        body_shape_type &body_shape = static_cast<body_shape_type&>(m_data.geometries.back());
                        detail::handle_body_geometry_config(geometry_tag,body_shape,stream);
                        m_data.surface.add_geometry(body_shape);
                    }
                    else if (geometry_type == "sine_shape")
                    {
                        m_data.geometries.push_back(new sine_shape_type());
                        sine_shape_type &sine_shape = static_cast<sine_shape_type&>(m_data.geometries.back());
                        detail::handle_sine_geometry_config(geometry_tag,sine_shape,stream);
                        m_data.surface.add_geometry(sine_shape);
                    }
                    else if (geometry_type == "circular_helix")
                    {
                        m_data.geometries.push_back(new circular_helix_type());
                        circular_helix_type &circular_helix = static_cast<circular_helix_type&>(m_data.geometries.back());
                        detail::handle_circular_helix_geometry_config(geometry_tag,circular_helix,stream);
                        m_data.surface.add_geometry(circular_helix);
                    }
                }

                return true;
            }

            bool warm_start(std::stringstream *stream)
            {
                fs::path data_path(m_data.data_path);

                if (fs::is_empty(data_path))
                    return false;

                if (stream)
                    *stream << "\t*************** Warm Starting the Simulation ***************" << std::endl;

                fs::directory_iterator it(data_path), end;

                size_t max = 0;

                while (it != end)
                {
                    size_t file_counter = get_file_number(fs::path(*it).filename());

                    if (file_counter >= max)
                    {
                        m_data.warm_starting_file = fs::path(*it).file_string();
                        max = file_counter;
                    }

                    ++it;
                }

                m_data.framecount = m_data.save_freq * max + 1;

                m_data.xmlcounter = max + 1;
                m_data.surface.time() = m_data.framecount *  m_data.timestep;

                if (stream)
                {
                    *stream << "\tSimiulation Time " << m_data.surface.time() << std::endl;
                    *stream << "\tSetting up frame " << m_data.xmlcounter << std::endl;
                    *stream << "\tProcessing timestep # " << m_data.framecount << std::endl;
                }

                if (stream)
                    *stream << "\tProcessing file " << m_data.warm_starting_file << std::endl;

                IO::read_vtk(0,m_data.surface,m_data.warm_starting_file);

                if (stream)
                    *stream << "\t**********************************************************\n" << std::endl;

//         if (path_array.size() > 1)
//         {
//             std::sort(path_array.begin(),path_array.end(),std::greater<fs::path>());
// //                     std::copy(path_array.begin(),path_array.end(), std::ostream_iterator<fs::path>(std::cout,"\n"));
// //                     std::cout << path_array[1].file_string() << std::endl;
//             m_data.warm_starting_file = path_array[1].file_string();
//             warm_start(m_data);
//         }
                return true;
            }

            size_t get_file_number(const std::string &in)
            {
                std::string out("");

                for (size_t i = 0; i < in.size(); ++i)
                    if (isdigit(in[i]))
                        out += in[i];

                return out.size()>0 ? boost::lexical_cast<size_t>(out) : 0;/*out;*/
            }


    };
}

#endif
// kate: indent-mode cstyle; space-indent on; indent-width 4;
