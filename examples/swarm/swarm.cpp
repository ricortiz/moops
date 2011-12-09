#include <iostream>
#include <iterator>
#include <memory>

#include "swarm.hpp"
#include "examples/type_binder.hpp"

class SwarmApp
#ifdef USE_QT_GUI
: public Gui<SwarmApp>
#endif
{
        enum
        {
            sdc_nodes = 5,
            M = 6,
            N = 100,
            head_particles = 20*12+1,
            num_particles = M * N + head_particles,
            data_size = 3 * num_particles,
            fmm_max_particles = 20,
            fmm_order = 3
        };
    protected:
        
#ifdef USE_CUDA_FLUID_SOLVER
        typedef float value_type;
#else
        typedef double value_type;
#endif
        typedef TypeBinder<value_type,num_particles> types;

    private:
        types::sine_type m_sine_geometry;
        types::swarm_surface_type m_swarm;
        types::swarm_boundary_type m_boundary;
        types::swarm_volume_type m_volume;
        value_type m_time_step;
        std::auto_ptr<types::vtk_writer> m_surface_writer;
        std::auto_ptr<types::vtk_writer> m_volume_writer;
        bool m_record;

    public:

        SwarmApp ( std::string &data_path ) :
                m_time_step ( .01 ),
                m_record ( true ),
                m_surface_writer ( new types::vtk_writer ( data_path + "/surface/" ) ),
                m_volume_writer ( new types::vtk_writer ( data_path + "/volume/" ) ),
                m_volume ( &m_boundary )
        {
            size_t *dims = m_sine_geometry.get_dimensions();
            value_type *x0 = m_sine_geometry.get_x0();
            dims[0] = M;
            dims[1] = N;
            dims[2] = 12;
            dims[3] = 21;
            x0[0] = x0[1] = x0[2] = 0;
            m_swarm.init_surface ( m_sine_geometry, m_boundary.positions(), num_particles );
//             m_swarm.init_volume ( m_sine_geometry, m_volume.positions() );
            m_swarm.set_springs ( m_boundary );
            m_boundary.set_surface ( m_swarm );
            m_boundary.fluid_solver().delta() = M_PI / N;

            m_sine_geometry.setCells();
            m_boundary.storage()->grid()->SetPolys ( m_sine_geometry.getCells() );

            m_surface_writer->setInput ( m_boundary.storage()->grid() );
//             m_volume_writer->setInput ( m_volume.storage()->grid() );
        }

        void run()
        {

            m_boundary.run ( m_time_step );
            m_volume.run ( m_time_step );
            if ( m_record )
            {
                m_surface_writer->write ( m_boundary.time() );
                m_volume_writer->write ( m_boundary.time() );
            }
        }

        void fake_run()
        {
//             particle_system_type::particle_type *particles = m_boundary.particles();
//             value_type dtheta = 2 * M_PI / M;
//             value_type dalpha = 2 * M_PI / N;
//             size_t lo, hi;
//             m_swarm.get_lohi ( lo, hi );
//             std::vector<value_type> &rscale = m_swarm.radius_scale ( m_boundary.time() );
//             for ( size_t i = lo, idx = lo * M; i < hi; ++i )
//                 for ( size_t j = 0; j < M; ++j, ++idx )
//                     m_sine_geometry.surface_point ( i, j, rscale[i-lo], particles[idx].position, dtheta, dalpha );
//             m_boundary.time() += m_time_step;
//             if ( m_record )
//                 m_surface_writer->write ( m_boundary.time() );
        }


    public:
        types::swarm_boundary_type *boundary()
        {
            return &m_boundary;
        }

};


int main ( int ac, char **av )
{
#ifdef USE_QT_GUI
    QVTKApplication app ( ac, av );
#endif

    std::string data_dir;
    if ( ac == 2 )
        data_dir = av[1];
    else
        return 0;

    SwarmApp swarm ( data_dir );
#ifdef USE_QT_GUI
    swarm.boundary()->storage()->grid()->PrintSelf ( std::cout, vtkIndent() );
    swarm.setActor ( swarm.boundary()->storage()->grid() );
    swarm.show();
    return app.exec();
#else
    while ( 1 )
        swarm.fake_run();
#endif
}
