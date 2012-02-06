#include <iostream>
#include <iterator>
#include <memory>

#include "examples/type_binder.hpp"
#include "valveless_heart.hpp"

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
        types::heart_pump_surface  m_heart_pump;
        value_type m_time_step;
        types::heart_vtk_storage m_vtk_storage;
        types::heart_vtk_writer m_surface_writer;
        bool m_record;

    public:

        ValvelessPump(std::string &data_path) :
                m_heart_pump(M, N),
                m_time_step(0.001),
                m_vtk_storage(m_heart_pump),
                m_surface_writer(data_path + "/surface/", m_vtk_storage, true),
                m_record(true)
        {
            m_heart_pump.fluid_solver().setDelta(.5*M_PI / N);
            m_heart_pump.geometry().setWaveSpeed(0.1);
            setCells();
        }

        void run()
        {
            if (m_record)
            {
                m_surface_writer.write(m_heart_pump.time());
//                 m_volume_writer->write(m_heart_pump.time());
            }
            m_heart_pump.run(m_time_step);
//             m_volume.run(m_time_step);
        }

        void fake_run()
        {
	    m_heart_pump.geometry().setWaveSpeed(200);
            if(m_record)
                m_surface_writer.write(m_heart_pump.time());
            types::particle_type *particles = m_heart_pump.particles();
	    value_type dtheta = 2 * M_PI / M;
            value_type dalpha = 2 * M_PI / N;
            size_t lo, hi;
            m_heart_pump.geometry().getForcingRange(lo, hi);
            const std::vector<value_type> &rscale = m_heart_pump.geometry().setPeristalticRadiusScaling(m_heart_pump.time());
            for (size_t i = lo, idx = lo * M; i < hi; ++i)
                for (size_t j = 0; j < M; ++j, ++idx)
                    m_heart_pump.geometry().surfacePoint(i, j, rscale[i - lo], particles[idx], dtheta, dalpha);
		
            m_heart_pump.time() += m_time_step;
        }
        void setCells()
        {
            m_vtk_storage.setInnerCells(M, N);
            m_vtk_storage.setEdgeCells(M, N);
            m_vtk_storage.setTopCells(M, N);
        }

        types::heart_vtk_storage &vtk_storage() { return m_vtk_storage; }
        types::heart_pump_surface &boundary() { return m_heart_pump; }
};


int main(int ac, char **av)
{
#ifdef USE_QT_GUI
    QVTKApplication app(ac, av);
#endif

    std::string data_dir;
    if (ac == 2)
        data_dir = av[1];
    else
        data_dir = "/home/rortiz/tmp/";

    ValvelessPump pump(data_dir);
#ifdef USE_QT_GUI
    pump.vtk_storage().grid()->PrintSelf(std::cout, vtkIndent());
    pump.setActor();
    pump.show();
    return app.exec();
#else
    int it = 0;
    double end_time = 1000;
    while (1)
    {
        pump.fake_run();
        logger.printTimer();
    }
    logger.writeTimer();
#endif
}
