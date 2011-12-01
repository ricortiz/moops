#include <iostream>
#include <iterator>
#include <memory>

#include <QVTKApplication.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkXMLPolyDataWriter.h>

#include "particle_system/time_integrator/sdc_integrator.hpp"
#include "particle_system/time_integrator/euler_integrator.hpp"
#include "particle_system/elastic_system/elastic_boundary.hpp"
#include "particle_system/elastic_system/spring_system.hpp"
#include "particle_system/particle_system.hpp"
#include "particle_system/fluid_solver/fmm_stokes_solver.hpp"
#ifdef CUDA_FLUID_SOLVER
#include "particle_system/fluid_solver/cuda_stokes_solver.hpp"
#endif
#include "particle_system/fluid_solver/direct_stokes_solver.hpp"
#include "math/ode_solver/sdc/explicit_sdc.hpp"
#include "math/ode_solver/euler/forward_euler.hpp"
#include "geometry/torus_geometry.hpp"
#include "geometry/oval_geometry.hpp"
#include "heart_pump.hpp"
#include "io/write_vtu.hpp"
#include "particle_system/particle_markers.hpp"
#include "gui/gui.hpp"

class ValvelessPump /*: public Gui<ValvelessPump>*/
{
    enum
    {
        sdc_nodes = 5,
        M = 20,
        N = 200,
        num_particles = M*N,
        data_size = 3*num_particles,
        fmm_max_particles = 20,
        fmm_order = 3
    };
protected:
#ifdef CUDA_FLUID_SOLVER
    typedef float value_type;
#else
    typedef double value_type;
#endif
    typedef ClenshawCurtis<value_type,data_size,sdc_nodes> integrator_type;
    typedef ExplicitSDC<value_type,integrator_type,sdc_nodes,sdc_nodes> sdc_type;
    typedef ForwardEuler<value_type> euler_type;
#ifdef CUDA_FLUID_SOLVER
    typedef vtkParticleSystemStorage<value_type,Particle<value_type>,PSYS::SURFACE,vtkFloatArray> surface_storage_type;
    typedef vtkParticleSystemStorage<value_type,Particle<value_type>,PSYS::VOLUME,vtkFloatArray> volume_storage_type;
    typedef ParticleSystem<value_type,sdc_type,num_particles,PSYS::SURFACE,surface_storage_type> particle_system_type;
    typedef ParticleSystem<value_type,euler_type,num_particles,PSYS::VOLUME,volume_storage_type> particle_system_tracers_type;
#else
    typedef ParticleSystem<value_type,sdc_type,num_particles> particle_system_type;
    typedef ParticleSystem<value_type,euler_type,num_particles,PSYS::VOLUME> particle_system_tracers_type;
#endif
    typedef HeartPump<value_type> heart_pump_surface_type;
    typedef TorusGeometry<value_type> torus_type;
    typedef OvalGeometry<value_type> oval_type;
    typedef SpringSystem<heart_pump_surface_type,particle_system_type> spring_system_type;
//         typedef FMMStokesSolver<particle_system_type,fmm_max_particles,fmm_order> fluid_solver_type;
#ifdef CUDA_FLUID_SOLVER
    typedef CudaStokesSolver<particle_system_type> fluid_solver_type;
#else
    typedef DirectStokesSolver<particle_system_type> fluid_solver_type;
#endif
    typedef ElasticBoundary<spring_system_type,fluid_solver_type,SDCIntegrator> boundary_type;
    typedef ParticleMarkers<particle_system_tracers_type,boundary_type,fluid_solver_type,EulerIntegrator> volume_type;

    typedef IO::VTKWriter<vtkXMLPolyDataWriter> vtk_writer;

private:
    oval_type m_oval_geometry;
    heart_pump_surface_type m_heart_pump;
    boundary_type m_boundary;
    volume_type m_volume;
    value_type m_time_step;
    std::auto_ptr<vtk_writer> m_surface_writer;
    std::auto_ptr<vtk_writer> m_volume_writer;
    bool m_record;

public:

    ValvelessPump(std::string &data_path) :
            m_time_step(.01),
            m_record(true),
            m_surface_writer(new vtk_writer(data_path+"/surface/")),
            m_volume_writer(new vtk_writer(data_path+"/volume/")),
            m_volume(&m_boundary)
    {
        size_t *dims = m_oval_geometry.get_dimensions();
        value_type *x0 = m_oval_geometry.get_x0();
        dims[0] = M;
        dims[1] = N;
        x0[0] = x0[1] = x0[2] = 0;
        m_heart_pump.init_surface(m_oval_geometry,m_boundary.positions());
        m_heart_pump.init_volume(m_oval_geometry,m_volume.positions());
        m_heart_pump.set_springs(m_boundary);
        m_boundary.set_surface(m_heart_pump);
        m_boundary.fluid_solver().delta() = M_PI/N;

        m_oval_geometry.setCells();
        m_boundary.storage()->grid()->SetPolys(m_oval_geometry.getCells());

        m_surface_writer->setInput(m_boundary.storage()->grid());
        m_volume_writer->setInput(m_volume.storage()->grid());
    }

    void run()
    {

        m_boundary.run(m_time_step);
        m_volume.run(m_time_step);
        if (m_record)
        {
            m_surface_writer->write(m_boundary.time());
            m_volume_writer->write(m_boundary.time());
        }
    }

    void fake_run()
    {
        particle_system_type::particle_type *particles = m_boundary.particles();
        size_t lo = 20, hi = 50;
        value_type dtheta = 2*M_PI/M;
        value_type dalpha = 2*M_PI/N;
        for(size_t i = lo, idx = lo*M; i < hi; ++i)
            for(size_t j = 0; j < M; ++j, ++idx)
                m_oval_geometry.surface_point(i,j,radius(m_boundary.time(),i-lo),particles[idx].position,dtheta,dalpha);
        m_boundary.time() += m_time_step;
        if (m_record)
        {
            m_surface_writer->write(m_boundary.time());
        }
    }

    value_type radius(value_type t, size_t i)
    {
        particle_system_type::particle_type *particles = m_boundary.particles();
        size_t lo = 20, hi = 50;
        value_type dtheta = 2*M_PI/M;
        value_type dalpha = 2*M_PI/N;
        
        // Range of stretchy part
        value_type x_1 = 0.0;
        value_type x_2 = 1.0;
        
        // Range of part to be squeezed
        value_type x_a = 0.2;
        value_type x_b = 0.8;
        
        // size of squeeze
        value_type Ls = 0.3;
        value_type r = .035;
        value_type deltar = m_oval_geometry.inner_radius() - r;
        
        value_type freq = 2.0;
        value_type wt = std::fmod(t/freq,1.0);
        value_type x_c = (1-wt)*x_a + wt*(x_b-Ls);
        value_type x_d = x_c + Ls;

        value_type dx = (x_2-x_1)/(hi-lo);

        std::vector<value_type> x(hi-lo), q(hi-lo);
        std::vector<value_type> filter(hi-lo);
        for(size_t k = 0; k < hi-lo; ++k)
        {
            x[k] = dx*k;            
            filter[k] = (x[k] > x_c) && (x[k] < x_d);
            q[k] = (x[k]-x_c)*(x[k]-x_d);
            q[k] *= q[k];
        }
        value_type scale = 0;
        for(size_t k = 0; k < hi-lo; ++k)
        {
            if(filter[k])
                if(q[k] > scale)
                    scale = q[k];
        }
        for(size_t k = 0; k < hi-lo; ++k)
            q[k] /= scale;

        return 1 - std::sin(M_PI*wt)*deltar*q[i]*filter[i]/m_oval_geometry.inner_radius();
    }

public:
    boundary_type *boundary()
    {
        return &m_boundary;
    }

};


int main(int ac, char **av)
{
//     QVTKApplication app(ac,av);

    std::string data_dir;
    if (ac == 2)
        data_dir = av[1];
    else
        return 0;

    ValvelessPump pump(data_dir);
//     pump.boundary()->storage()->grid()->PrintSelf(std::cout, vtkIndent());
//     pump.setActor(pump.boundary()->storage()->grid());
//     pump.show();
    while (1)
        pump.run();
//     return app.exec();
}
