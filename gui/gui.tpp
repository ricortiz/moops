/*=========================================================================
 *
 *  Program:   Visualization Toolkit
 *  Module:    GUI4.cxx
 *
 *  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 *  All rights reserved.
 *  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.
 *
 *     This software is distributed WITHOUT ANY WARRANTY; without even
 *     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 *     PURPOSE.  See the above copyright notice for more information.
 *
 * =========================================================================*/

/*=========================================================================
 *
 *  Copyright 2004 Sandia Corporation.
 *  Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 *  license for use of this work by or on behalf of the
 *  U.S. Government. Redistribution and use in source and binary forms, with
 *  or without modification, are permitted provided that this Notice and any
 *  statement of authorship are reproduced on all copies.
 *
 * =========================================================================*/

/*========================================================================
 * For general information about using VTK and Qt, see:
 * http://www.trolltech.com/products/3rdparty/vtksupport.html
 * =========================================================================*/

#include "gui.hpp"

#include <QMenu>
#include <QLabel>
#include <QApplication>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenuBar>
#include <QStatusBar>

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCommand.h>
#include <vtkEventQtSlotConnect.h>
#include <vtkPolyDataMapper.h>
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

template<typename app>
Gui<app>::Gui() : GuiBase()
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
    
    const double angleSensitivity=0.01;
    const double translationSensitivity=0.001;
    
    QVTKInteractor *interactor = m_vtk_widget->GetInteractor();
    vtkInteractorStyle *s = static_cast<vtkInteractorStyle *>(interactor->GetInteractorStyle());
    vtkTDxInteractorStyleCamera *t= static_cast<vtkTDxInteractorStyleCamera *>(s->GetTDxStyle());
    
    t->GetSettings()->SetAngleSensitivity(angleSensitivity);
    t->GetSettings()->SetTranslationXSensitivity(translationSensitivity);
    t->GetSettings()->SetTranslationYSensitivity(translationSensitivity);
    t->GetSettings()->SetTranslationZSensitivity(translationSensitivity);
    
    // add a renderer
    m_vtk_renderer = vtkRenderer::New();
    m_vtk_renderer->SetBackground(84./110,89./110,109./110);
    m_vtk_widget->GetRenderWindow()->AddRenderer(m_vtk_renderer);
    
    // add a popup menu for the window and connect it to our slot
    QMenu* popupMenu = new QMenu(m_vtk_widget);
    popupMenu->addAction("Reset");
    connect(popupMenu, SIGNAL(triggered(QAction*)), this, SLOT(reset(QAction*)));

    m_connections = vtkEventQtSlotConnect::New();    
    // get right mouse pressed with high priority
    m_connections->Connect(m_vtk_widget->GetRenderWindow()->GetInteractor(),vtkCommand::RightButtonPressEvent,this,SLOT(popup(vtkObject*, unsigned long, void*, void*, vtkCommand*)),popupMenu, 1.0);
    // update coords as we move through the window
    m_connections->Connect(m_vtk_widget->GetRenderWindow()->GetInteractor(),vtkCommand::MouseMoveEvent,this,SLOT(updateCoords(vtkObject*)));    
    m_connections->PrintSelf(cout, vtkIndent());
}

template<typename app>
Gui<app>::~Gui()
{
    m_vtk_renderer->Delete();
    m_connections->Delete();
}

template<typename app>
void Gui<app>::popup(vtkObject * obj, unsigned long,void * client_data, void *,vtkCommand * command)
{
    // get interactor
    vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::SafeDownCast(obj);
    // consume event so the interactor style doesn't get it
    command->AbortFlagOn();
    // get popup menu
    QMenu* popupMenu = static_cast<QMenu*>(client_data);
    // get event location
    int* sz = interactor->GetSize();
    int* position = interactor->GetEventPosition();
    // remember to flip y
    QPoint pt = QPoint(position[0], sz[1]-position[1]);
    // map to global
    QPoint global_pt = popupMenu->parentWidget()->mapToGlobal(pt);
    // show popup menu at global point
    popupMenu->popup(global_pt);
}

template<typename app>
void Gui<app>::setActor(vtkPolyData *poly_data)
{
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoother = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoother->SetInput(poly_data);
    smoother->SetNumberOfIterations(1);
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputConnection(smoother->GetOutputPort());
    normals->SetSplitting(1);
    normals->SetNonManifoldTraversal(1);
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    mapper->SetInputConnection(normals->GetOutputPort());
    actor->SetMapper(mapper);
    m_vtk_renderer->AddViewProp(actor);
}

template<typename app>
void Gui<app>::updateCoords(vtkObject* obj)
{
    // get interactor
    vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::SafeDownCast(obj);
    // get event position
    int event_pos[2];
    interactor->GetEventPosition(event_pos);
    // update label
    QString str;
    str.sprintf("x=%d : y=%d", event_pos[0], event_pos[1]);
    coord->setText(str);
}
