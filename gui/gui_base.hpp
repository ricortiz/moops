#ifndef GUI_BASE_HPP
#define GUI_BASE_HPP
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

#include <QMainWindow>

class QAction;
class QVBoxLayout;
class QHBoxLayout;
class QLabel;
class QMenuBar;
class QMenu;
class QStatusBar;
class QTimer;
class QRadialGradient;
class QMouseEvent;
class QTextDocument;
class QPainter;
class vtkObject;
class vtkCommand;

class GuiBase : public QMainWindow
{
        Q_OBJECT

    public:
        GuiBase();

    public slots:
        void popup(vtkObject * obj, unsigned long, void * client_data, void *, vtkCommand * command);
        void updateCoords(vtkObject* obj);
        void reset(QAction*);

    protected:
        QVBoxLayout     *vboxLayout;
        QHBoxLayout     *hboxLayout;
        QWidget         *centralwidget;
        QAction         *actionExit;
        QLabel          *coord;
        QMenuBar        *menubar;
        QMenu           *menuFile;
        QStatusBar      *statusbar;
        QTimer          *m_timer;
};

#endif
