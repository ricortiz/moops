#include <iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cstdlib>
#ifdef USE_QT_GUI
#include <QApplication>
#endif
#include "math/fluid_solver/stokes/fmm/hybrid/tree_renderer.hpp"

#include "utils/logger.hpp"

Logger logger;
#include "math/fluid_solver/stokes/hybrid_fmm_stokes_solver.hpp"

template<typename value_type>
struct random_generator
{
    value_type operator()()
    {
        return rand()/((value_type)RAND_MAX+1);
    }
};

int hybrid_fmm_stokes_solver(int ac,char **av)
{
    srand(0);
    typedef float value_type;
    size_t num_sources = 1 << 10;
    size_t num_targets = 1 << 10;
    size_t size_sources = 3*num_sources;
    size_t size_targets = 3*num_targets;
    std::vector<value_type> sources(size_sources), targets(size_targets), velocities(size_targets), forces(size_sources);
    value_type delta = .01;
    value_type domain[2][3] = {{0,0,0},{10,10,10}};
    HybridFmmStokesSolver<value_type> fmm(num_sources);
    std::generate(sources.begin(),sources.end(),random_generator<value_type>());
    std::generate(targets.begin(),targets.end(),random_generator<value_type>());
    std::generate(forces.begin(),forces.end(),random_generator<value_type>());    
    std::cout << "sources = ["; std::copy(sources.begin(),sources.end(),std::ostream_iterator<value_type>(std::cout," ")); std::cout << "];" << std::endl;
    std::cout << "targets = ["; std::copy(targets.begin(),targets.end(),std::ostream_iterator<value_type>(std::cout," ")); std::cout << "];" << std::endl;
    std::cout << "forces = ["; std::copy(forces.begin(),forces.end(),std::ostream_iterator<value_type>(std::cout," ")); std::cout << "];" << std::endl;
    fmm.setDelta(delta);
    fmm(0,&sources[0],&velocities[0],&forces[0]);

    std::fill(velocities.begin(),velocities.end(),0.0);
    fmm.allPairs();
#ifdef USE_QT_GUI
    QApplication app(ac, av);

    if (!QGLFormat::hasOpenGL())
    {
        std::cerr << "This system has no OpenGL support" << std::endl;
        return 1;
    }
    QGL::setPreferredPaintEngine(QPaintEngine::OpenGL);
    OctreeRenderer tree_renderer;
    tree_renderer.init(octree);
    tree_renderer.setWindowTitle(QObject::tr("Quad Tree"));
    tree_renderer.setMinimumSize(200, 200);
    tree_renderer.resize(800, 600);
    tree_renderer.show();
    return app.exec();
#else
    return 0;    
#endif
}