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

#include "gui_base.hpp"

#include <QMenu>
#include <QTimer>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenuBar>
#include <QStatusBar>
#include <QApplication>
#include <QPainter>

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

void GuiBase::popup(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command)
{
    vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::SafeDownCast(obj); // get interactor
    command->AbortFlagOn();          // consume event so the interactor style doesn't get it
    QMenu* popupMenu = static_cast<QMenu*>(client_data);       // get popup menu
    int* sz = interactor->GetSize();
    int* position = interactor->GetEventPosition();       // get event location
    QPoint pt = QPoint(position[0], sz[1] - position[1]);     // remember to flip y
    QPoint global_pt = popupMenu->parentWidget()->mapToGlobal(pt);     // map to global
    popupMenu->popup(global_pt);         // show popup menu at global point
}

void GuiBase::updateCoords(vtkObject* obj)
{
    vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::SafeDownCast(obj); // get interactor
    int event_pos[2];
    interactor->GetEventPosition(event_pos);        // get event position
    QString str;
    str.sprintf("x=%d : y=%d", event_pos[0], event_pos[1]);
    coord->setText(str);          // update label
}


