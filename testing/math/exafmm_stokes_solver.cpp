#include <iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cstdlib>
#ifdef USE_QT_GUI
#include <QApplication>
#include "math/fluid_solver/stokes/fmm/hybrid/tree_renderer.hpp"
#endif

#include "utils/vtk_storage_wrapper.hpp"
#include "io/write_vtu.hpp"

#include "math/fluid_solver/stokes/exafmm_stokes_solver.hpp"

template<typename value_type>
struct random_generator
{
    value_type operator()()
    {
        return rand() / ((value_type)RAND_MAX + 1);
    }
};

int exafmm_stokes_solver(int ac, char **av)
{
    srand(0);
    bool dump_tree_data = false;
    vtkSmartPointer<vtkPolyData>   poly_data = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints>     points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkFloatArray> data_points = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkCellArray>  cells = vtkSmartPointer<vtkCellArray>::New();
    typedef float value_type;
    size_t num_sources = 1 << 8;
    size_t num_targets = 1 << 8;
    size_t size_sources = 3 * num_sources;
    size_t size_targets = 3 * num_targets;
    std::vector<value_type> sources(size_sources), targets(size_targets), velocities(size_targets), forces(size_sources);
    value_type delta = .1;
    ExaFmmStokesSolver<value_type> fmm(num_sources);
    std::generate(sources.begin(), sources.end(), random_generator<value_type>());
    std::generate(targets.begin(), targets.end(), random_generator<value_type>());
    std::generate(forces.begin(), forces.end(), random_generator<value_type>());
    std::cout << "sources = ["; std::copy(sources.begin(), sources.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "targets = ["; std::copy(targets.begin(), targets.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "forces = ["; std::copy(forces.begin(), forces.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    fmm.setDelta(delta);
    fmm(0, &sources[0], &velocities[0], &forces[0]);
    std::cout << "velocities = ["; std::copy(velocities.begin(), velocities.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
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
