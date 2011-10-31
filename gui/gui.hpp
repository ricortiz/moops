#ifndef GUI_HPP
#define GUI_HPP

#include <QMainWindow>
#include <vtkSmartPointer.h>

class vtkPolyData;
class vtkRenderer;
class vtkEventQtSlotConnect;
class vtkObject;
class vtkCommand;
class QVTKWidget;

class QAction;
class QVBoxLayout;
class QHBoxLayout;
class QLabel;
class QMenuBar;
class QMenu;
class QStatusBar;


class GUI : public QMainWindow
{
    Q_OBJECT
public:
    GUI();
    ~GUI();

    void SetActor(vtkPolyData *);
    
public slots:
    void updateCoords(vtkObject*);
    void popup(vtkObject * obj, unsigned long,void * client_data, void *,vtkCommand * command);
    void reset(QAction*);
    void setup();
    void retranslateUi();


protected:
    vtkRenderer* m_vtk_renderer;
    vtkEventQtSlotConnect* Connections;

    QWidget *centralwidget;
    QAction *actionExit;
    QAction *actionE_xit;
    QVBoxLayout *vboxLayout;
    QHBoxLayout *hboxLayout;
    QVTKWidget *m_vtk_widget;
    QLabel *coord;
    QMenuBar *menubar;
    QMenu *menuFile;
    QMenu *menuFile_2;
    QMenu *menu_File;
    QStatusBar *statusbar;

};

#endif
