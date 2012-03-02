#ifndef GUI_HPP
#define GUI_HPP
/****************************************************************************
** MOOPS -- Modular Object Oriented Particle Simulator
** Copyright (C) 2011-2012  Ricardo Ortiz <ortiz@unc.edu>
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
****************************************************************************/

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
