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


class GuiBase : public QMainWindow
{
        Q_OBJECT
        
    public:
        GuiBase();


    public slots:
        void reset(QAction*);


    protected:

        QWidget *centralwidget;
        QAction *actionExit;
        QAction *actionE_xit;
        QVBoxLayout *vboxLayout;
        QHBoxLayout *hboxLayout;
        QLabel *coord;
        QMenuBar *menubar;
        QMenu *menuFile;
        QMenu *menuFile_2;
        QMenu *menu_File;
        QStatusBar *statusbar;

};

#endif
