#ifndef GUI_BASE_HPP
#define GUI_BASE_HPP

#include <QMainWindow>

class QAction;
class QVBoxLayout;
class QHBoxLayout;
class QLabel;
class QMenuBar;
class QMenu;
class QStatusBar;
class vtkObject;
class vtkCommand;

class GuiBase : public QMainWindow
{
        Q_OBJECT
        
    public:
        GuiBase();


    public slots:
    void popup(vtkObject * obj, unsigned long,void * client_data, void *,vtkCommand * command);
    void updateCoords(vtkObject* obj);
        void reset(QAction*);


    protected:

        QWidget *centralwidget;
        QAction *actionExit;
        QVBoxLayout *vboxLayout;
        QHBoxLayout *hboxLayout;
        QLabel *coord;
        QMenuBar *menubar;
        QMenu *menuFile;
        QStatusBar *statusbar;

};

#endif
