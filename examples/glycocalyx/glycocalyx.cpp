#include <iostream>
#include <iterator>
#include <memory>

#include "utils/type_binder.hpp"

class GlycocalyxApp
#if defined USE_QT_GUI
            : public Gui<GlycocalyxApp>
#elif defined USE_PV_COPROCESSOR
            : public ParaviewCoprocessor<GlycocalyxApp>
#endif
{
    public:
        std::string name() { return "Swarm"; }

    protected:
        enum
        {
            num_towers = 49,
            M = 6,
            N = 20,
            tower_particles = M*N,
            num_particles = num_towers * tower_particles,
            data_size = 3 * num_particles
        };
        typedef float value_type;
        typedef TypeBinder<value_type> types;

        types::glycocalyx_surface m_glycocalyx;
        types::glycocalyx_vtk_storage m_vtk_storage;
        types::vtk_poly_writer m_surface_writer;
        bool m_record;

    public:

        GlycocalyxApp() :
                m_glycocalyx(M, N, num_towers), m_vtk_storage(m_glycocalyx), m_record(true)
        {
            m_glycocalyx.fluid_solver().setDelta(.02);
            for (int i = 0; i < num_towers; ++i)
                setCells(i*tower_particles);
        }

        int init(int ac, char **av)
        {
            if (ac < 2)
            {
                std::cout << "Usage: " << av[0] << " <data path>" << std::endl;
                return 0;
            }
            std::string path = av[1];
            m_surface_writer.setDataPath(path + "surface/glycocalyx");
            return 1;
        }

        void run()
        {
            if (m_record)
                write();
            m_glycocalyx.run();
        }

        void setCells(size_t offset = 0)
        {
            m_vtk_storage.setInnerCells(M, N - 1, offset);
            m_vtk_storage.setEdgeCells(M, N - 1, offset);
        }

        types::glycocalyx_vtk_storage &vtk_storage() { return m_vtk_storage; }

        void write()
        {
            m_surface_writer.write(m_vtk_storage.grid(), m_glycocalyx.time());
        }
    public:
        types::glycocalyx_surface &boundary()
        {
            return m_glycocalyx;
        }
};


int main(int ac, char **av)
{
    GlycocalyxApp glycocalyx;
    glycocalyx.init(ac, av);
#if defined USE_QT_GUI
    return runGui(ac, av, glycocalyx);
#elif defined USE_PV_COPROCESSOR
    glycocalyx.initCoprocessor(ac, av);
    return runCoprocessor(glycocalyx);
#else
    while (1)
        glycocalyx.run();
#endif
    return 0;
}

