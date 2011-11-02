#ifndef GUI_HPP
#define GUI_HPP

#include "gui_base.hpp"

class vtkPolyData;
class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;
class QVTKWidget;
class GuiBase;

template<typename app>
class Gui : public GuiBase
{
public:
    Gui();
    ~Gui();
    void setActor(vtkPolyData *);
public slots:
    void popup(vtkObject * obj, unsigned long,void * client_data, void *,vtkCommand * command);
    void updateCoords(vtkObject* obj);
    
protected:
    app *derived() {
        return static_cast<app*>(this);
    }

private:

    vtkRenderer* m_vtk_renderer;
    vtkEventQtSlotConnect* m_connections;
    QVTKWidget *m_vtk_widget;
};

#include "gui.tpp"

#endif
