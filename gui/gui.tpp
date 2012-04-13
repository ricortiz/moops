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

#include <QMenu>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenuBar>
#include <QStatusBar>
#include <QGLWidget>
#include <QRadialGradient>
#include <QMouseEvent>
#include <QTextDocument>

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCommand.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkActor.h>
#include <vtkInteractorStyle.h>
#include <vtkTDxInteractorStyleCamera.h>
#include <vtkTDxInteractorStyleSettings.h>
#include <vtkSmartPointer.h>
#include <vtkDecimatePro.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPolyDataNormals.h>
#include <QVTKInteractor.h>
#include <QVTKWidget.h>

template<typename app_type>
Gui<app_type>::Gui() : GuiBase()
{
    m_vtk_widget = new QVTKWidget(centralwidget);
    m_vtk_widget->setObjectName(QString::fromUtf8("m_vtk_widget"));
    hboxLayout->addWidget(m_vtk_widget);

    vtkSmartPointer<vtkRenderWindow> renwin = vtkSmartPointer<vtkRenderWindow>::New();

    // Activate 3DConnexion device
#ifdef VTK_USE_TDX
    m_vtk_widget->SetUseTDx(true);
#endif

    m_vtk_widget->SetRenderWindow(renwin);

    const double angleSensitivity = 0.1;
    const double translationSensitivity = 0.01;

    QVTKInteractor *interactor = m_vtk_widget->GetInteractor();
    vtkInteractorStyle *s = static_cast<vtkInteractorStyle *>(interactor->GetInteractorStyle());
    vtkTDxInteractorStyleCamera *t = static_cast<vtkTDxInteractorStyleCamera *>(s->GetTDxStyle());

    t->GetSettings()->SetAngleSensitivity(angleSensitivity);
    t->GetSettings()->SetTranslationXSensitivity(translationSensitivity);
    t->GetSettings()->SetTranslationYSensitivity(translationSensitivity);
    t->GetSettings()->SetTranslationZSensitivity(translationSensitivity);

    // add a renderer
    m_vtk_renderer = vtkRenderer::New();
    m_vtk_renderer->SetBackground(105. / 255, 146. / 255, 182. / 255);
    m_vtk_renderer->SetBackground2(16. / 255, 56. / 255, 121. / 255);
    m_vtk_renderer->SetGradientBackground(true);
    m_vtk_widget->GetRenderWindow()->AddRenderer(m_vtk_renderer);

    // add a popup menu for the window and connect it to our slot
    QMenu* popupMenu = new QMenu(m_vtk_widget);
    popupMenu->addAction("Reset");
    connect(popupMenu, SIGNAL(triggered(QAction*)), this, SLOT(reset(QAction*)));

    m_connections = vtkEventQtSlotConnect::New();
    // get right mouse pressed with high priority
    m_connections->Connect(m_vtk_widget->GetRenderWindow()->GetInteractor(), vtkCommand::RightButtonPressEvent, this, SLOT(popup(vtkObject*, unsigned long, void*, void*, vtkCommand*)), popupMenu, 1.0);
    // update coords as we move through the window
    m_connections->Connect(m_vtk_widget->GetRenderWindow()->GetInteractor(), vtkCommand::MouseMoveEvent, this, SLOT(updateCoords(vtkObject*)));
//     m_connections->PrintSelf(cout, vtkIndent());
}

template<typename app_type>
Gui<app_type>::~Gui()
{
    m_vtk_renderer->Delete();
    m_connections->Delete();
}


template<typename app_type>
void Gui<app_type>::setGridActor(bool smooth = false)
{
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    if (smooth)
    {
        vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        smoother->SetInput(app()->vtk_storage().grid());
        smoother->SetNumberOfIterations(15);
        vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        normals->SetInputConnection(smoother->GetOutputPort());
        normals->SetSplitting(1);
        normals->SetNonManifoldTraversal(1);
        mapper->SetInputConnection(normals->GetOutputPort());
    }
    else
    {
        mapper->SetInput(app()->vtk_storage().grid());
    }
    actor->SetMapper(mapper);

    m_vtk_renderer->AddViewProp(actor);
}

template<typename app_type>
void Gui<app_type>::setBoxActor()
{
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    mapper->SetInput(app()->vtk_storage().box());
    actor->SetMapper(mapper);
    actor->GetProperty()->SetOpacity(.1);
    actor->GetProperty()->SetEdgeVisibility(1);

    m_vtk_renderer->AddViewProp(actor);
}

