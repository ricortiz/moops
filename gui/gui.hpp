#ifndef GUI_HPP
#define GUI_HPP

#include <QVTKApplication.h>
#include "gui_base.hpp"

class vtkPolyData;
class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;
class QVTKWidget;

template<typename app_type>
class Gui : public GuiBase
{
    public:
        Gui();
        ~Gui();
        void setGridActor(bool);
        void setBoxActor();
    public slots:
        void popup(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command);
        void updateCoords(vtkObject* obj);

    protected:
        app_type *app()
        {
            return static_cast<app_type*>(this);
        }

    private:

        vtkRenderer* m_vtk_renderer;
        vtkEventQtSlotConnect* m_connections;
        QVTKWidget *m_vtk_widget;
};

#include "gui.tpp"

template<typename app_type>
int runGui(Gui<app_type> &app)
{
    QVTKApplication app_window(ac, av);
    app.setGridActor(true);
    app.setBoxActor();
    app.show();
    return app_window.exec();
}

#endif
