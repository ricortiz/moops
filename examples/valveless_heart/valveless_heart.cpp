#include <iostream>
#include <iterator>
#include <memory>

#include "utils/logger.hpp"

Logger logger;

#include "valveless_heart.hpp"
#include "examples/type_binder.hpp"

class ValvelessPump
#ifdef USE_QT_GUI
    : public Gui<ValvelessPump>
#endif
{
        enum
        {
            sdc_nodes = 5,
            M = 20,
            N = 200,
            num_particles = M * N,
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
        typedef TypeBinder<value_type> types;

    private:
        types::oval_type m_oval_geometry;
        types::heart_pump_surface_type m_heart_pump;
        types::heart_pump_boundary_type m_boundary;
        types::heart_pump_volume_type m_volume;
        value_type m_time_step;
        std::auto_ptr<types::vtk_writer> m_surface_writer;
        std::auto_ptr<types::vtk_writer> m_volume_writer;
        bool m_record;

    public:

        ValvelessPump(std::string &data_path) : m_boundary(data_size),
            m_time_step(.01),
            m_record(true),
            m_surface_writer(new types::vtk_writer(data_path + "/surface/")),
            m_volume_writer(new types::vtk_writer(data_path + "/volume/")),
            m_volume(m_boundary, data_size)
        {
            m_oval_geometry.setDimensions(M,N);
            m_oval_geometry.setX0(0,0,0);
            m_oval_geometry.setForcingRange(10,50);
            m_oval_geometry.setWaveSpeed(.001);
            m_heart_pump.initSurface(m_oval_geometry, m_boundary.particles());
//         std::cout << "p = [";
//         std::copy(m_boundary.positions(), m_boundary.positions() + data_size, std::ostream_iterator<value_type>(std::cout, ","));
//         std::cout << "];\n";
            m_heart_pump.initVolume(m_oval_geometry, m_volume.particles());
            m_heart_pump.setSprings(m_boundary);
            m_boundary.setSurface(m_heart_pump);
            m_boundary.ode_rhs().delta() = .5*M_PI/N;

            m_oval_geometry.setCells();
            m_boundary.storage()->grid()->SetPolys(m_oval_geometry.getCells());

            m_surface_writer->setInput(m_boundary.storage()->grid());
            m_volume_writer->setInput(m_volume.storage()->grid());
        }

        void run()
        {
            if(m_record)
            {
                m_surface_writer->write(m_boundary.time());
                m_volume_writer->write(m_boundary.time());
            }
            m_boundary.run(m_time_step);
            m_volume.run(m_time_step);
        }

        void fake_run()
        {
            types::particle_system_type::particle_type *particles = m_boundary.particles();
            value_type dtheta = 2 * M_PI / M;
            value_type dalpha = 2 * M_PI / N;
            size_t lo, hi;
            m_oval_geometry.getForcingRange(lo, hi);
            const std::vector<value_type> &rscale = m_oval_geometry.setRadiusScaling(m_boundary.time());
// 	    std::copy(rscale.begin(),rscale.end(),std::ostream_iterator<value_type>(std::cout, " ")); std::cout << std::endl;
            for(size_t i = lo, idx = lo * M; i < hi; ++i)
                for(size_t j = 0; j < M; ++j, ++idx)
                    m_oval_geometry.surface_point(i, j, rscale[i - lo], particles[idx], dtheta, dalpha);
            m_boundary.time() += m_time_step;
            if(m_record)
                m_surface_writer->write(m_boundary.time());
        }


    public:
        types::heart_pump_boundary_type *boundary()
        {
            return &m_boundary;
        }

        value_type timestep() { return m_time_step; }

};


int main(int ac, char **av)
{
#ifdef USE_QT_GUI
    QVTKApplication app(ac, av);
#endif

    std::string data_dir;
    if(ac == 2)
        data_dir = av[1];
    else
        return 0;

    ValvelessPump pump(data_dir);
#ifdef USE_QT_GUI
    pump.boundary()->storage()->grid()->PrintSelf(std::cout, vtkIndent());
    pump.setActor(pump.boundary()->storage()->grid());
    pump.show();
    return app.exec();
#else
    int it = 0;
    double end_time = 1000;
    while(it++*pump.timestep() < end_time)
    {
        pump.run();
        logger.printTimer();
    }
    logger.writeTimer();
#endif
}
