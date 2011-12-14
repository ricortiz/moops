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

#include "gui_base.hpp"

#include <QMenu>
#include <QTimer>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenuBar>
#include <QStatusBar>
#include <QApplication>

#include <vtkRenderWindowInteractor.h>
#include <vtkCommand.h>

GuiBase::GuiBase()
{
    this->setWindowTitle(QApplication::translate("GuiBase", "MainWindow", 0, QApplication::UnicodeUTF8));

    if (this->objectName().isEmpty())
        this->setObjectName(QString::fromUtf8("Gui"));
    this->resize(800, 600);
    
    actionExit    = new QAction(this);
    centralwidget = new QWidget(this);
    vboxLayout    = new QVBoxLayout(centralwidget);
    hboxLayout    = new QHBoxLayout();
    coord         = new QLabel(centralwidget);
    menubar       = new QMenuBar(this);
    menuFile      = new QMenu(menubar);
    statusbar     = new QStatusBar(this);
    m_timer       = new QTimer(this);
    
    actionExit->setObjectName(QString::fromUtf8("actionExit"));
    centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
    vboxLayout->setObjectName(QString::fromUtf8("vboxLayout"));
    hboxLayout->setObjectName(QString::fromUtf8("hboxLayout"));
    coord->setObjectName(QString::fromUtf8("coord"));
    menubar->setObjectName(QString::fromUtf8("menubar"));
    menuFile->setObjectName(QString::fromUtf8("menuFile"));
    statusbar->setObjectName(QString::fromUtf8("statusbar"));
    
    vboxLayout->addLayout(hboxLayout);
    coord->setAlignment(Qt::AlignCenter);
    vboxLayout->addWidget(coord);
    menubar->setGeometry(QRect(0, 0, 442, 29));
    menubar->addAction(menuFile->menuAction());
    menuFile->addAction(actionExit);
    actionExit->setText(QApplication::translate("GuiBase", "E&xit", 0, QApplication::UnicodeUTF8));
    menuFile->setTitle(QApplication::translate("GuiBase", "&File", 0, QApplication::UnicodeUTF8));
    coord->setText(QApplication::translate("GuiBase", "TextLabel", 0, QApplication::UnicodeUTF8));
    this->setCentralWidget(centralwidget);
    this->setMenuBar(menubar);
    this->setStatusBar(statusbar);
    
    QObject::connect(actionExit, SIGNAL(triggered()), this, SLOT(close()));
    QMetaObject::connectSlotsByName(this);
}


void GuiBase::reset(QAction* color)
{
    //     m_vtk_renderer->ResetCamera();
    //     m_vtk_widget->update();
}

void GuiBase::popup(vtkObject * obj, unsigned long,void * client_data, void *,vtkCommand * command)
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
