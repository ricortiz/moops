#include <iostream>
#include <iterator>
#include <memory>

#include "examples/type_binder.hpp"

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
            num_sperms = 25,
            Mt = 6,
            Nt = 100,
            Mh = 12,
            Nh = 21,
            tail_particles = Mt * Nt,
            head_particles = Mh * (Nh - 1) + 1,
            total_particles = tail_particles + head_particles,
            num_particles = num_sperms * total_particles,
            data_size = 3 * num_particles
        };
        typedef float value_type;
        typedef TypeBinder<value_type> types;

        types::glycocalyx_surface m_glycocalyx;
        types::value_type m_time_step;
        types::glycocalyx_vtk_storage m_vtk_storage;
        types::vtk_poly_writer m_surface_writer;
        types::vtk_unstructured_writer m_octree_writer;
        bool m_record;

    public:

        GlycocalyxApp() :
                m_glycocalyx(Mt, Nt, Mh, Nh, num_sperms),
                m_time_step(0.01),
                m_vtk_storage(m_glycocalyx),
                m_record(true)
        {
            m_glycocalyx.fluid_solver().setDelta(.02);
            m_glycocalyx.geometry().setWaveSpeed(0.01);
            for (int i = 0; i < num_sperms; ++i)
                setCells(i*total_particles);
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
            m_octree_writer.setDataPath(path + "octree/glycocalyx");
            return 1;
        }

        void run()
        {
            if (m_record)
                write();
            m_glycocalyx.run(m_time_step);
        }

        void fake_run()
        {
            m_glycocalyx.geometry().setWaveSpeed(20);
            if (m_record)
                write();
            types::particle_type *particles = m_glycocalyx.particles();
            value_type dtheta = 2 * M_PI / Mt;
            value_type dalpha = 2 * M_PI / Nt;
            value_type x0[3];
            m_glycocalyx.geometry().getX0(x0);
            for (int s = 0; s < num_sperms; ++s)
                for (int i = 0, idx = head_particles + s * total_particles; i < Nt; ++i)
                    for (int j = 0; j < Mt; ++j, ++idx)
                    {
                        m_glycocalyx.geometry().surfacePoint(i, j, m_glycocalyx.time(), particles[idx], dtheta, dalpha);
                        //                         particles[idx].position[0] += x0[0];
                        particles[idx].position[1] += x0[1] + (1.0 - 1.0 / 42);
                        //                         particles[idx].position[2] += x0[2];
                    }
            m_glycocalyx.time() += m_time_step;
        }

        void setCells(size_t offset = 0)
        {
            // Set cells for the nose of the head
            {
                for (size_t i = 1; i < Mh; ++i)
                {
                    vtkIdType cell[3] = {0 + offset, i + 1 + offset, i + offset};
                    m_vtk_storage.cells()->InsertNextCell(3, cell);
                }
                vtkIdType cell[3] = {0 + offset, 1 + offset, Mh + offset};
                m_vtk_storage.cells()->InsertNextCell(3, cell);
            }
            m_vtk_storage.setInnerCells(Mh, Nh - 1, offset + 1);
            m_vtk_storage.setEdgeCells(Mh, Nh - 1, offset + 1);
            setJunctionCells(offset);
            m_vtk_storage.setInnerCells(Mt, Nt, offset + head_particles);
            m_vtk_storage.setEdgeCells(Mt, Nt, offset + head_particles);
        }

        void setJunctionCells(size_t offset)
        {
            int factor = Mh / Mt;
            int head_points = head_particles - Mh;
            for (size_t i = 0; i < Mt; ++i)
            {
                vtkIdType cells[3][3] = {{i + head_particles + offset, head_points + i *factor + 1 + offset, head_points + i *factor + offset},
                    {i + head_particles + offset, (i + 1) % Mt + head_particles + offset, head_points + i *factor + 1 + offset},
                    {(i + 1) % Mt + head_particles + offset, head_points + (i *factor + 2) % Mh + offset, head_points + i *factor + 1 + offset}
                };
                for (int k = 0; k < 3; ++k)
                    m_vtk_storage.cells()->InsertNextCell(3, cells[k]);
            }
        }

        types::glycocalyx_vtk_storage &vtk_storage() { return m_vtk_storage; }

        void write()
        {
            m_surface_writer.write(m_vtk_storage.grid(), m_glycocalyx.time());
            vtkOctreeStorage vtk_octree;
            vtk_octree.setOctree(m_glycocalyx.fluid_solver());
            m_octree_writer.write(vtk_octree.getBox(), m_glycocalyx.time());
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

