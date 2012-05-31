
#include "utils/type_binder.hpp"
#include <mgl2/mgl.h>


class App
{
protected:
    enum
    {
        M = 20,
        N = 200,
        num_sperms = 1,
        Mt = 6,
        Nt = 100,
        Mh = 12,
        Nh = 21,
    };
    typedef double value_type;
    typedef TypeBinder<value_type> types;

    types::swarm_surface            m_swarm;
    types::heart_pump_surface       m_heart;

public:
    App() : m_swarm ( Mt, Nt, Mh, Nh, num_sperms ), m_heart ( M, N )
    {
        m_heart.fluid_solver().setDelta ( .005 );
        m_heart.geometry().setWaveSpeed ( .01 );
        m_swarm.fluid_solver().setDelta ( .02 );
        m_swarm.geometry().setWaveSpeed ( .01 );
    }

    types::heart_pump_surface &HeartSolver()
    {
        return m_heart;
    }

    types::swarm_surface &SwarmSolver()
    {
        return m_swarm;
    }


};

template<typename value_type>
void PlotConvergence ( mglGraph &graph, std::vector<value_type> &x_axis, std::vector<value_type> &y_axis )
{
    graph.Title ( "Convergence plot and rate" );
    graph.SetOrigin(0,0);
    mglData x ,y;
    x.Set(&x_axis[0],x_axis.size());
    y.Set(&y_axis[0],y_axis.size());
//     y.Linear ( float ( y_axis.size()-1 ) );
//     x.Linear ( float ( x_axis.size()-1 ) );
    graph.Axis(); 
    graph.Label('x',"dt");
    graph.Label('y',"error");
    graph.SetRanges(x.Minimal()-.02,x.Maximal()+.02,y.Minimal()-.02,y.Maximal()+.02);
    graph.SetFunc ( "lg(x)","lg(y)" );
    graph.Plot ( x, y, "or-" ); 
    graph.WriteFrame("convergence.eps");
}


template<typename stepper_type, typename value_type>
void RunStepper ( stepper_type &solver, int iterations, value_type total_time )
{
    value_type dt = total_time / iterations;
    std::cout << "Computing solution with dt = " << dt << "; num_iterations = " << iterations << std::endl;
    for ( int i = 0; i < iterations; ++i )
        solver.run ( dt );
}

template<typename stepper_type>
int ComputeConvergenceRate ( stepper_type &solver, std::vector<int> &num_iterations )
{
    typedef typename stepper_type::value_type value_type;
    std::vector<value_type> solution0 ( solver.data_size() ), reference_solution ( solver.data_size() );
    std::copy ( solver.positions(),solver.positions() +solver.data_size(),solution0.begin() );
    const value_type total_time = .1;
    const int reference_iterations = 64;
    RunStepper ( solver, reference_iterations, total_time );
    std::copy ( solver.positions(),solver.positions() +solver.data_size(),reference_solution.begin() );
    solver.clearTime();
    solver.clearForces();
    solver.clearVelocities();
    std::copy ( solution0.begin(),solution0.end(),solver.positions() );

    std::vector<value_type> dt ( num_iterations.size(),total_time );
    std::vector<int>::const_iterator ni = num_iterations.begin();
    typename std::vector<value_type>::iterator i = dt.begin();
    for ( ; i != dt.end(); ++i,++ni )
        *i /= *ni;

    std::vector<value_type> e ( num_iterations.size(),0.0 );
    for ( int i = 0; i < num_iterations.size(); ++i )
    {
        RunStepper ( solver, num_iterations[i], total_time );

        std::vector<value_type> error ( solver.data_size(), 0.0 );
        std::transform ( reference_solution.begin(), reference_solution.end(), solver.positions(), error.begin(), std::minus<value_type>() );
        e[i] = std::sqrt(std::inner_product ( error.begin(), error.end(), error.begin(), 0.0 ));
        solver.clearTime();
        solver.clearForces();
        solver.clearVelocities();
        std::copy ( solution0.begin(),solution0.end(),solver.positions() );
    }
    mglGraph graph;
    PlotConvergence ( graph,dt,e );
}


int main()
{
//     App app;
//     std::vector<int> num_iterations;
//     num_iterations.push_back ( 21 );
//     num_iterations.push_back ( 13 );
//     num_iterations.push_back ( 8 );
//     num_iterations.push_back ( 5 );
// //     ComputeConvergenceRate(app.HeartSolver(),num_iterations);
//     ComputeConvergenceRate ( app.SwarmSolver(),num_iterations );
    mglGraph graph;
    std::vector<double> dt,e;
    dt.push_back(.0001);
    dt.push_back(.001);
    dt.push_back(.01);
    dt.push_back(.1);
    e.push_back(.0001);
    e.push_back(.001);
    e.push_back(.01);
    e.push_back(.1);
    PlotConvergence ( graph,dt,e );
}
