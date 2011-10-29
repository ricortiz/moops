#ifndef PARTICLE_RENDERER_HPP
#define PARTICLE_RENDERER_HPP

#include <cmath>
#include <QtOpenGL/QGLWidget>
#include <QtGui/QRadialGradient>
#include <QtGui/QMouseEvent>
#include <QtGui/QTextDocument>


class VowelCube : public QGLWidget
{
        Q_OBJECT

    public:
        VowelCube(QWidget *parent = 0): QGLWidget(parent)
        {
            setFormat(QGLFormat(QGL::SampleBuffers));

            m_rotationX = -38.0;
            m_rotationY = -58.0;
            m_rotationZ = 0.0;
            m_scaling = 1.0;

            
        }
        ~VowelCube()
        {
            makeCurrent();
            glDeleteLists(m_glObject, 1);
        }
        template<typename tree_type>
        void init( tree_type *tree )
        {
            createGradient();
            renderBox(tree->root());
        }
    protected:
        void paintEvent(QPaintEvent * /* event */)
        {
            QPainter painter(this);
            drawBackground(&painter);
            drawCube();
//             drawLegend(&painter);
        }
        void mousePressEvent(QMouseEvent *event)
        {
            m_lastPos = event->pos();
        }
        void mouseMoveEvent(QMouseEvent *event)
        {
            GLfloat dx = GLfloat(event->x() - m_lastPos.x()) / width();
            GLfloat dy = GLfloat(event->y() - m_lastPos.y()) / height();

            if (event->buttons() & Qt::LeftButton)
            {
                m_rotationX += 180 * dy;
                m_rotationY += 180 * dx;
                update();
            }
            else if (event->buttons() & Qt::RightButton)
            {
                m_rotationX += 180 * dy;
                m_rotationZ += 180 * dx;
                update();
            }
            m_lastPos = event->pos();
        }
        void wheelEvent(QWheelEvent *event)
        {
            double numDegrees = -event->delta() / 8.0;
            double numSteps = numDegrees / 15.0;
            m_scaling *= std::pow(1.125, numSteps);
            update();
        }

    private:
        void createGradient()
        {
            m_gradient.setCoordinateMode(QGradient::ObjectBoundingMode);
            m_gradient.setCenter(0.45, 0.50);
            m_gradient.setFocalPoint(0.40, 0.45);
            m_gradient.setColorAt(0.0, QColor(105, 146, 182));
            m_gradient.setColorAt(0.4, QColor(81, 113, 150));
            m_gradient.setColorAt(0.8, QColor(16, 56, 121));
        }

        template<typename box_type>
        void renderBox(const box_type &box)
        {
            makeCurrent();

            glShadeModel(GL_FLAT);

            m_glObject = glGenLists(1);
            glNewList(m_glObject, GL_COMPILE);
            qglColor(QColor(255, 239, 191));
            glLineWidth(1.0);

            createBoxes(box);

            glEndList();
        }

        template<typename box_type>
        void createBoxes(box_type *box)
        {
            typedef typename box_type::box_iterator box_iterator;
            typedef typename box_type::partilce_iterator particle_iterator;

            glBegin(GL_LINES);
            glVertex3f(box->x()+box->u(), box->y()+box->v(), box->z()-box->w());
            glVertex3f(box->x()-box->u(), box->y()+box->v(), box->z()-box->w());
            glVertex3f(box->x()+box->u(), box->y()-box->v(), box->z()-box->w());
            glVertex3f(box->x()-box->u(), box->y()-box->v(), box->z()-box->w());
            glVertex3f(box->x()+box->u(), box->y()-box->v(), box->z()+box->w());
            glVertex3f(box->x()-box->u(), box->y()-box->v(), box->z()+box->w());
            glEnd();

            glBegin(GL_LINE_LOOP);
            glVertex3f(box->x()+box->u(), box->y()+box->v(), box->z()+box->w());
            glVertex3f(box->x()+box->u(), box->y()+box->v(), box->z()-box->w());
            glVertex3f(box->x()+box->u(), box->y()-box->v(), box->z()-box->w());
            glVertex3f(box->x()+box->u(), box->y()-box->v(), box->z()+box->w());
            glVertex3f(box->x()+box->u(), box->y()+box->v(), box->z()+box->w());
            glVertex3f(box->x()-box->u(), box->y()+box->v(), box->z()+box->w());
            glVertex3f(box->x()-box->u(), box->y()+box->v(), box->z()-box->w());
            glVertex3f(box->x()-box->u(), box->y()-box->v(), box->z()-box->w());
            glVertex3f(box->x()-box->u(), box->y()-box->v(), box->z()+box->w());
            glVertex3f(box->x()-box->u(), box->y()+box->v(), box->z()+box->w());
            glEnd();

            if (box->children().size() == 0)
                for (particle_iterator p = box->particles().begin(), end = box->particles().end(); p != end; ++p)
                {
                    GLUquadric * qobj = gluNewQuadric();
                    glPushMatrix();
                    glTranslatef(*p->position[0], *p->position[1], *p->position[2]);
                    GLint slices = 8;
                    GLint stacks = 8;
                    gluSphere(qobj, .01, slices, stacks);
                    glPopMatrix();
                    gluDeleteQuadric(qobj);
                }

            for (box_iterator b = box->children().begin(), end = box->children().end(); b != end; ++b)
            {
#pragma omp task firstprivate(b)
                createBoxes(*b);
            }
#pragma omp taskwait

        }

        void drawBackground(QPainter *painter)
        {
            painter->setPen(Qt::NoPen);
            painter->setBrush(m_gradient);
            painter->drawRect(rect());
        }
        void drawCube()
        {
            glPushAttrib(GL_ALL_ATTRIB_BITS);

            glMatrixMode(GL_PROJECTION);
            glPushMatrix();
            glLoadIdentity();
            GLfloat x = 3.0 * GLfloat(width()) / height();
            glOrtho(-x, +x, -3.0, +3.0, 4.0, 15.0);

            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glTranslatef(0.0, 0.0, -10.0);
            glScalef(m_scaling, m_scaling, m_scaling);

            glRotatef(m_rotationX, 1.0, 0.0, 0.0);
            glRotatef(m_rotationY, 0.0, 1.0, 0.0);
            glRotatef(m_rotationZ, 0.0, 0.0, 1.0);

            glEnable(GL_MULTISAMPLE);

            glCallList(m_glObject);

//             setFont(QFont("Times", 24));
//             qglColor(QColor(255, 223, 127));
// 
//             renderText(+1.1, +1.1, +1.1, QChar('a'));
//             renderText(-1.1, +1.1, +1.1, QChar('e'));
//             renderText(+1.1, +1.1, -1.1, QChar('o'));
//             renderText(-1.1, +1.1, -1.1, QChar(0x00F6));
//             renderText(+1.1, -1.1, +1.1, QChar(0x0131));
//             renderText(-1.1, -1.1, +1.1, QChar('i'));
//             renderText(+1.1, -1.1, -1.1, QChar('u'));
//             renderText(-1.1, -1.1, -1.1, QChar(0x00FC));

            glMatrixMode(GL_MODELVIEW);
            glPopMatrix();

            glMatrixMode(GL_PROJECTION);
            glPopMatrix();

            glPopAttrib();
        }
        void drawLegend(QPainter *painter)
        {
            const int Margin = 11;
            const int Padding = 6;

            QTextDocument textDocument;
            textDocument.setDefaultStyleSheet("* { color: #FFEFEF }");
            textDocument.setHtml("<h4 align=\"center\">Vowel Categories</h4>"
                                 "<p align=\"center\"><table width=\"100%\">"
                                 "<tr><td>Open:<td>a<td>e<td>o<td>&ouml;"
                                 "<tr><td>Close:<td>&#305;<td>i<td>u<td>&uuml;"
                                 "<tr><td>Front:<td>e<td>i<td>&ouml;<td>&uuml;"
                                 "<tr><td>Back:<td>a<td>&#305;<td>o<td>u"
                                 "<tr><td>Round:<td>o<td>&ouml;<td>u<td>&uuml;"
                                 "<tr><td>Unround:<td>a<td>e<td>&#305;<td>i"
                                 "</table>");
            textDocument.setTextWidth(textDocument.size().width());

            QRect rect(QPoint(0, 0), textDocument.size().toSize()
                       + QSize(2 * Padding, 2 * Padding));
            painter->translate(width() - rect.width() - Margin,
                               height() - rect.height() - Margin);
            painter->setPen(QColor(255, 239, 239));
            painter->setBrush(QColor(255, 0, 0, 31));
            painter->drawRect(rect);

            painter->translate(Padding, Padding);
            textDocument.drawContents(painter);
        }

        GLuint m_glObject;
        QRadialGradient m_gradient;
        GLfloat m_rotationX;
        GLfloat m_rotationY;
        GLfloat m_rotationZ;
        GLfloat m_scaling;
        QPoint m_lastPos;
};



#endif  // PARTICLE_RENDERER_HPP

