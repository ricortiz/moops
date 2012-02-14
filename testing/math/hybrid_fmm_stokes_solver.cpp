#include <iostream>
#include<vector>
#include<algorithm>
#include<iterator>
#include<cstdlib>
#ifdef USE_QT_GUI
#include <QApplication>
#include "math/fluid_solver/stokes/fmm/hybrid/tree_renderer.hpp"
#endif

#include "utils/logger.hpp"
#include "utils/vtk_storage_wrapper.hpp"
#include "io/write_vtu.hpp"

Logger logger;
#include "math/fluid_solver/stokes/hybrid_fmm_stokes_solver.hpp"

template<typename value_type>
struct random_generator
{
    value_type operator()()
    {
        return rand() / ((value_type)RAND_MAX + 1);
    }
};
void write_vtk_node(Node &node, const std::string &name)
{

    vtk_octree_storage vtk_octree;
    IO::VtkWriter<vtkXMLUnstructuredGridWriter> writer(name);
    float points[8][3] = {{0}};
    vtk_octree.getPoints(&node, octree.edge_length, points);
    vtk_octree.addBox(points);
    vtk_octree.setCurrentPoints();
    writer.write(vtk_octree.getBox(), 0, false);

}

void write_vtk_node_points(Node &node, const std::string &name)
{
    vtkSmartPointer<vtkPolyData>        poly_data = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkFloatArray>      data_points = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkPoints>          points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray>       cells = vtkSmartPointer<vtkCellArray>::New();
    IO::VtkWriter<vtkXMLPolyDataWriter> writer(name);
    // TODO THis is wrong
    data_points->SetArray(octree.bodies[node.pArrayLow].position, 3*(node.pArrayHigh - node.pArrayLow + 1), 1);
    data_points->SetNumberOfComponents(3);
    data_points->SetName("positions");

    points->SetData(data_points);
    poly_data->SetPoints(points);
    for (vtkIdType i = 0; i <= node.pArrayHigh - node.pArrayLow; ++i)
        cells->InsertNextCell(VTK_VERTEX, &i);
    poly_data->SetVerts(cells);
    writer.write(poly_data, 0, false);
}

void write_vtk_octree()
{
    for (int i = 0; i < octree.rootInfo->subtree_size; ++i)
    {
//         if(i == 331)
        if (!octree.nodes[i].isParent && octree.nodes[i].pArrayLow >= 0 && octree.nodes[i].pArrayHigh >= 0)
        {
            std::stringstream s;
            s << "data/node_" << i;
            write_vtk_node(octree.nodes[i], s.str());
            for (int j = 0; j < *octree.nodes[i].list1Count; ++j)
            {
                std::stringstream n, p;
                n << "data/list1_" << "node_" << i << "_id_" << j;
                write_vtk_node(*octree.nodes[i].list1[j], n.str());
                n << "_points";
//                 write_vtk_node_points(*octree.nodes[i].list1[j],n.str());
            }
            for (int j = 0; j < *octree.nodes[i].list2Count; ++j)
            {
                std::stringstream n;
                n << "data/list2_" << "node_" << i << "_id_" << j;
                write_vtk_node(*octree.nodes[i].list2[j], n.str());
                n << "_points";
//                 write_vtk_node_points(*octree.nodes[i].list2[j],n.str());
            }
            for (int j = 0; j < *octree.nodes[i].list3Count; ++j)
            {
                std::stringstream n;
                n << "data/list3_" << "node_" << i << "_id_" << j;
                write_vtk_node(*octree.nodes[i].list3[j], n.str());
                n << "_points";
//                 write_vtk_node_points(*octree.nodes[i].list3[j],n.str());
            }
            for (int j = 0; j < *octree.nodes[i].list4Count; ++j)
            {
                std::stringstream n;
                n << "data/list4_" << "node_" << i << "_id_" << j;
                write_vtk_node(*octree.nodes[i].list4[j], n.str());
                n << "_points";
//                 write_vtk_node_points(*octree.nodes[i].list4[j],n.str());
            }
            for (int j = 0; j < *octree.nodes[i].colleagueCount; ++j)
            {
                std::stringstream n;
                n << "data/colleagues_" << "node_" << i << "_id_" << j;
                write_vtk_node(*octree.nodes[i].Colleagues[j], n.str());
                n << "_points";
//                 write_vtk_node_points(*octree.nodes[i].Colleagues[j],n.str());
            }
        }
    }

}

int hybrid_fmm_stokes_solver(int ac, char **av)
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
    value_type delta = .0001;
    HybridFmmStokesSolver<value_type> fmm(num_sources);
    std::generate(sources.begin(), sources.end(), random_generator<value_type>());
    std::generate(targets.begin(), targets.end(), random_generator<value_type>());
    std::generate(forces.begin(), forces.end(), random_generator<value_type>());
    std::cout << "sources = ["; std::copy(sources.begin(), sources.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "targets = ["; std::copy(targets.begin(), targets.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    std::cout << "forces = ["; std::copy(forces.begin(), forces.end(), std::ostream_iterator<value_type>(std::cout, " ")); std::cout << "];" << std::endl;
    fmm.setDelta(delta);
    fmm(0, &sources[0], &velocities[0], &forces[0]);

    std::fill(velocities.begin(), velocities.end(), 0.0);
    fmm.allPairs();
    if (dump_tree_data)
    {
        IO::VtkWriter<vtkXMLPolyDataWriter> writer("data/points");
        data_points->SetArray(&sources[0], sources.size(), 1);
        data_points->SetNumberOfComponents(3);
        data_points->SetName("positions");
        points->SetData(data_points);
        poly_data->SetPoints(points);
        for (vtkIdType i = 0; i < num_sources; ++i)
            cells->InsertNextCell(VTK_VERTEX, &i);
        poly_data->SetVerts(cells);
        writer.write(poly_data, 0, false);
        write_vtk_octree();
    }
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
