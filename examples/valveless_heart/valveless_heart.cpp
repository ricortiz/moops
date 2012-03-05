#include <iostream>
#include <iterator>
#include <memory>

#include "examples/type_binder.hpp"
#include "valveless_heart.hpp"

class ValvelessPump
#if defined USE_QT_GUI
            : public Gui<ValvelessPump>
#elif defined USE_PV_COPROCESSOR
            : public ParaviewCoprocessor<ValvelessPump>
#endif
{
    public:
        std::string name() { return "Heart Pump"; }

    protected:
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
        typedef float value_type;
        typedef TypeBinder<value_type> types;

        types::heart_pump_surface       m_heart_pump;
        value_type                      m_time_step;
        types::heart_vtk_storage        m_vtk_storage;
        types::vtk_poly_writer          m_surface_writer;
        bool                            m_record;

    public:

        ValvelessPump() :
                m_heart_pump(M, N),
                m_time_step(0.01),
                m_vtk_storage(m_heart_pump),
                m_record(true)
        {
            m_heart_pump.fluid_solver().setDelta(.5*M_PI / N);
            m_heart_pump.geometry().setWaveSpeed(m_time_step+.5*m_time_step);
            setCells();
        }

        int init(int ac, char **av)
        {
            if (ac < 2)
            {
                std::cout << "Usage: " << av[0] << " <data path>" << std::endl;
                return 0;
            }
            std::string path = av[1];
            m_surface_writer.setDataPath(path + "/surface/heart");
            return 1;
        }
        
        void run()
        {
            if (m_record)
                write();
            m_heart_pump.run(m_time_step);
        }

        void fake_run()
        {
            m_heart_pump.geometry().setWaveSpeed(200);
            if (m_record)
                write();
            types::particle_type *particles = m_heart_pump.particles();
            value_type dtheta = 2 * M_PI / M;
            value_type dalpha = 2 * M_PI / N;
            m_heart_pump.computeForces(m_heart_pump.time());
            std::vector<value_type> rscale;
            m_heart_pump.geometry().getRadiusScaling(rscale);
            size_t lo, hi;
            m_heart_pump.geometry().getForcingRange(lo, hi);
            for (size_t i = lo, idx = lo * M; i < hi; ++i)
                for (size_t j = 0; j < M; ++j, ++idx)
                {
                    m_heart_pump.geometry().surfacePoint(i, j, rscale[i - lo], particles[idx], dtheta, dalpha);
                    particles[idx].position[2] += .15;
                }
                
            m_heart_pump.time() += m_time_step;
        }
        
        void setCells()
        {
            m_vtk_storage.setInnerCells(M, N);
            m_vtk_storage.setEdgeCells(M, N);
            m_vtk_storage.setTopCells(M, N);
        }
        
        types::heart_vtk_storage &vtk_storage() { return m_vtk_storage; }
        void write(){ m_surface_writer.write(m_vtk_storage.grid(), m_heart_pump.time()); }
        types::heart_pump_surface &boundary() { return m_heart_pump; }
};


int main(int ac, char **av)
{
#if defined USE_QT_GUI
    QVTKApplication app(ac, av);
    ValvelessPump heart;
    if(!heart.init(ac, av))
        return 1;
    heart.setGridActor(true);
    heart.setBoxActor();
    heart.show();
    return app.exec();
#elif defined USE_PV_COPROCESSOR
    ValvelessPump heart;
    if(!heart.init(ac, av))
        return 1;
    heart.initCoprocessor(ac, av);
    return runCoprocessor(heart);
#else
    ValvelessPump heart;
    if(!heart.init(ac, av))
        return 1;
    while (1)
        heart.run();
    return 0;
#endif
}
