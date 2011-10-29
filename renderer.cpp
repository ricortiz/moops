#include <QtGui/QApplication>
#include <iostream>

#include "io/particle_renderer.hpp"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    
    if (!QGLFormat::hasOpenGL()) {
        std::cerr << "This system has no OpenGL support" << std::endl;
        return 1;
    }
    QGL::setPreferredPaintEngine(QPaintEngine::OpenGL);
    VowelCube cube;
    cube.setWindowTitle(QObject::tr("Vowel Cube"));
    cube.setMinimumSize(200, 200);
    cube.resize(450, 450);
    cube.show();
    
    return app.exec();
}