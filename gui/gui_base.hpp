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

        QVBoxLayout  *vboxLayout;
        QHBoxLayout  *hboxLayout;
        QWidget *centralwidget;
        QAction *actionExit;
        QLabel  *coord;
        QMenuBar *menubar;
        QMenu   *menuFile;
        QStatusBar *statusbar;
        QTimer *m_timer;
};

#endif
