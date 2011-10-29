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

GuiBase::GuiBase()
{
    this->setup();

    // create a window to make it stereo capable and give it to QVTKWidget
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

    // put cone in one window
    Connections = vtkEventQtSlotConnect::New();

    // get right mouse pressed with high priority
    Connections->Connect(m_vtk_widget->GetRenderWindow()->GetInteractor(),vtkCommand::RightButtonPressEvent,this,SLOT(popup(vtkObject*, unsigned long, void*, void*, vtkCommand*)),popupMenu, 1.0);

    // update coords as we move through the window
    Connections->Connect(m_vtk_widget->GetRenderWindow()->GetInteractor(),vtkCommand::MouseMoveEvent,this,SLOT(updateCoords(vtkObject*)));

    Connections->PrintSelf(cout, vtkIndent());

}

void GuiBase::SetActor(vtkPolyData *poly_data)
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

GuiBase::~GuiBase()
{
    m_vtk_renderer->Delete();
    Connections->Delete();
}


void GuiBase::updateCoords(vtkObject* obj)
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

void GuiBase::popup(vtkObject * obj, unsigned long,void * client_data, void *,vtkCommand * command)
{
    // A note about context menus in Qt and the QVTKWidget
    // You may find it easy to just do context menus on right button up,
    // due to the event proxy mechanism in place.

    // That usually works, except in some cases.
    // One case is where you capture context menu events that
    // child windows don't process.  You could end up with a second
    // context menu after the first one.

    // See QVTKWidget::ContextMenuEvent enum which was added after the
    // writing of this example.

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

void GuiBase::reset(QAction* color)
{
//     m_vtk_renderer->ResetCamera();
//     m_vtk_widget->update();
}


void GuiBase::retranslateUi()
{
    this->setWindowTitle(QApplication::translate("GuiBase", "MainWindow", 0, QApplication::UnicodeUTF8));
    actionExit->setText(QApplication::translate("GuiBase", "Exit", 0, QApplication::UnicodeUTF8));
    actionE_xit->setText(QApplication::translate("GuiBase", "E&xit", 0, QApplication::UnicodeUTF8));
    coord->setText(QApplication::translate("GuiBase", "TextLabel", 0, QApplication::UnicodeUTF8));
    menuFile->setTitle(QApplication::translate("GuiBase", "File", 0, QApplication::UnicodeUTF8));
    menuFile_2->setTitle(QApplication::translate("GuiBase", "File", 0, QApplication::UnicodeUTF8));
    menu_File->setTitle(QApplication::translate("GuiBase", "&File", 0, QApplication::UnicodeUTF8));
}

void GuiBase::setup()
{
    if (this->objectName().isEmpty())
        this->setObjectName(QString::fromUtf8("GuiBase"));
    this->resize(442, 361);
    actionExit = new QAction(this);
    actionExit->setObjectName(QString::fromUtf8("actionExit"));
    actionE_xit = new QAction(this);
    actionE_xit->setObjectName(QString::fromUtf8("actionE_xit"));
    centralwidget = new QWidget(this);
    centralwidget->setObjectName(QString::fromUtf8("centralwidget"));

    vboxLayout = new QVBoxLayout(centralwidget);
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    hboxLayout = new QHBoxLayout();
    hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));

    m_vtk_widget = new QVTKWidget(centralwidget);
    m_vtk_widget->setObjectName(QString::fromUtf8("m_vtk_widget"));
    QSizePolicy sizePolicy(static_cast<QSizePolicy::Policy>(7), static_cast<QSizePolicy::Policy>(3));
    sizePolicy.setHorizontalStretch(0);
    sizePolicy.setVerticalStretch(0);
    sizePolicy.setHeightForWidth(m_vtk_widget->sizePolicy().hasHeightForWidth());
    m_vtk_widget->setSizePolicy(sizePolicy);

    hboxLayout->addWidget(m_vtk_widget);


    vboxLayout->addLayout(hboxLayout);

    coord = new QLabel(centralwidget);
    coord->setObjectName(QString::fromUtf8("coord"));
    coord->setAlignment(Qt::AlignCenter);

    vboxLayout->addWidget(coord);

    this->setCentralWidget(centralwidget);
    menubar = new QMenuBar(this);
    menubar->setObjectName(QString::fromUtf8("menubar"));
    menubar->setGeometry(QRect(0, 0, 442, 29));
    menuFile = new QMenu(menubar);
    menuFile->setObjectName(QString::fromUtf8("menuFile"));
    menuFile_2 = new QMenu(menubar);
    menuFile_2->setObjectName(QString::fromUtf8("menuFile_2"));
    menu_File = new QMenu(menubar);
    menu_File->setObjectName(QString::fromUtf8("menu_File"));
    this->setMenuBar(menubar);
    statusbar = new QStatusBar(this);
    statusbar->setObjectName(QString::fromUtf8("statusbar"));
    this->setStatusBar(statusbar);

    menubar->addAction(menu_File->menuAction());
    menu_File->addAction(actionE_xit);

//     retranslateUi();
    QObject::connect(actionE_xit, SIGNAL(triggered()), this, SLOT(close()));

    QMetaObject::connectSlotsByName(this);
}

