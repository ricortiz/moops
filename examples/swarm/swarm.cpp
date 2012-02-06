#include <iostream>
#include <iterator>
#include <memory>

#include "examples/type_binder.hpp"

class SwarmApp
#ifdef USE_QT_GUI
    : public Gui<SwarmApp>
#endif
{
        enum
        {
            num_sperms = 1,
            Mt = 6,
            Nt = 100,
            Mh = 12,
            Nh = 21,
            tail_particles = Mt*Nt,
            head_particles = Mh*(Nh-1) + 1,
            total_particles = tail_particles + head_particles,
            num_particles = num_sperms * total_particles,
            data_size = 3 * num_particles
        };
    protected:
        typedef double value_type;
        typedef TypeBinder<value_type> types;
        
    private:
        types::swarm_surface m_swarm;
        types::value_type m_time_step;
        types::swarm_vtk_storage m_vtk_storage;
        types::swarm_vtk_writer m_surface_writer;
        bool m_record;

    public:

        SwarmApp(std::string &data_path) :
            m_swarm(Mt,Nt,Mh,Nh,num_sperms),
            m_time_step(0.001),
            m_vtk_storage(m_swarm),
            m_surface_writer(data_path + "/surface/",m_vtk_storage,true),
            m_record(true)
            
        {
            m_swarm.fluid_solver().setDelta(.02);
            m_swarm.geometry().setWaveSpeed(0.1);
            for(int i = 0; i < num_sperms; ++i)
                setCells(i*total_particles);
        }

        void run()
        {
            if(m_record)
                m_surface_writer.write(m_swarm.time());
            m_swarm.run(m_time_step);
        }

        void fake_run()
        {
            m_swarm.geometry().setWaveSpeed(20);
            if(m_record)
                m_surface_writer.write(m_swarm.time());
            types::particle_type *particles = m_swarm.particles();
            value_type dtheta = 2 * M_PI / Mt;
            value_type dalpha = 2 * M_PI / Nt;
            value_type x0[3];
            m_swarm.geometry().getX0(x0);
            for(int s = 0; s < num_sperms; ++s)
                for(int i = 0, idx = head_particles + s * total_particles; i < Nt; ++i)
                    for(int j = 0; j < Mt; ++j, ++idx)
                    {
                        m_swarm.geometry().surfacePoint(i, j, m_swarm.time(), particles[idx], dtheta, dalpha);
                        particles[idx].position[0] += x0[0];
                        particles[idx].position[1] += x0[1] + (1.0 - 1.0 / 42);
                        particles[idx].position[2] += x0[2];
                    }
            m_swarm.time() += m_time_step;
        }

        void setCells(size_t offset = 0)
        {
            // Set cells for the nose of the head
            {
                for(size_t i = 1; i < Mh; ++i)
                {
                    vtkIdType cell[3] = {0 + offset, i + 1 + offset, i + offset};
                    m_vtk_storage.cells()->InsertNextCell(3, cell);
                }
                vtkIdType cell[3] = {0 + offset, 1 + offset, Mh + offset};
                m_vtk_storage.cells()->InsertNextCell(3, cell);
            }
            m_vtk_storage.setInnerCells(Mh,Nh-1,offset+1);
            m_vtk_storage.setEdgeCells(Mh,Nh-1,offset+1);
            setJunctionCells(offset);
            m_vtk_storage.setInnerCells(Mt,Nt,offset+head_particles);
            m_vtk_storage.setEdgeCells(Mt,Nt,offset+head_particles);
        }

        void setJunctionCells(size_t offset)
        {
            int factor = Mh / Mt;
            int head_points = head_particles - Mh;
            for(size_t i = 0; i < Mt; ++i)
            {
                vtkIdType cells[3][3] = {{i + head_particles + offset, head_points + i *factor + 1 + offset, head_points + i *factor + offset},
                    {i + head_particles + offset, (i + 1) % Mt + head_particles + offset, head_points + i *factor + 1 + offset},
                    {(i + 1) % Mt + head_particles + offset, head_points + (i *factor + 2) % Mh + offset, head_points + i *factor + 1 + offset}
                };
                for(int k = 0; k < 3; ++k)
                    m_vtk_storage.cells()->InsertNextCell(3, cells[k]);
            }
        }

        types::swarm_vtk_storage &vtk_storage() { return m_vtk_storage; }
        
    public:
        types::swarm_surface &boundary()
        {
            return m_swarm;
        }

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

    SwarmApp swarm(data_dir);
#ifdef USE_QT_GUI
    swarm.vtk_storage().grid()->PrintSelf(std::cout, vtkIndent());
    swarm.setActor(swarm.vtk_storage().grid());
    swarm.show();
    return app.exec();
#else
    while(1)
        swarm.run();
#endif
}
