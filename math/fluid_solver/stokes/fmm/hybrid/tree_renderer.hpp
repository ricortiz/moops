#ifndef OCTREE_RENDERER_HPP
#define OCTREE_RENDERER_HPP

#include <cmath>
#include <QtOpenGL/QGLWidget>
#include <QtGui/QRadialGradient>
#include <QtGui/QMouseEvent>
#include <QtGui/QTextDocument>
#include <iostream>

class OctreeRenderer : public QGLWidget
{
        Q_OBJECT

    public:
        OctreeRenderer(QWidget *parent = 0): QGLWidget(parent), num_particles(0), num_boxes(0)
        {
            setFormat(QGLFormat(QGL::SampleBuffers));

            m_rotationX = -38.0;
            m_rotationY = -58.0;
            m_rotationZ = 0.0;
            m_scaling = .3;


        }
        ~OctreeRenderer()
        {
            makeCurrent();
            glDeleteLists(m_glObject, 1);
        }
        template<typename box_array_type>
        void init(box_array_type &boxes)
        {
            createGradient();
            renderBox(boxes.root,boxes.bodies,boxes.edge_length[0]);
            num_boxes = 0;//boxes.size();
//             depth = tree->depth();
//             num_particles = 0;//boxes[0].particles().size();
//             max_particles = tree->get_max_particles();
        }
    protected:
        void paintEvent(QPaintEvent * /* event */)
        {
            QPainter painter(this);
            drawBackground(&painter);
            drawCube();
            drawLegend(&painter);
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
            m_gradient.setCenter(0.0, 0.0);
            m_gradient.setFocalPoint(0.01, 0.01);
            m_gradient.setColorAt(0.0, QColor(105, 146, 182));
            m_gradient.setColorAt(0.4, QColor(81, 113, 150));
            m_gradient.setColorAt(0.8, QColor(16, 56, 121));
        }

        template<typename box_type, typename particle_array>
        void renderBox(const box_type *box, particle_array *particles, double extent)
        {
            makeCurrent();

            glShadeModel(GL_FLAT);

            m_glObject = glGenLists(1);
            glNewList(m_glObject, GL_COMPILE);
            qglColor(QColor(255, 239, 191));
            glLineWidth(1.0);

            createBoxes(box,particles,extent);

            glEndList();
        }

        template<typename box_type, typename particle_array>
        void createBoxes(box_type *box, particle_array* particles, double extent)
        {
            {
                glBegin(GL_LINES);
                glVertex3f(box->mid_x+extent, box->mid_y+extent, box->mid_z-extent);
                glVertex3f(box->mid_x-extent, box->mid_y+extent, box->mid_z-extent);
                glVertex3f(box->mid_x+extent, box->mid_y-extent, box->mid_z-extent);
                glVertex3f(box->mid_x-extent, box->mid_y-extent, box->mid_z-extent);
                glVertex3f(box->mid_x+extent, box->mid_y-extent, box->mid_z+extent);
                glVertex3f(box->mid_x-extent, box->mid_y-extent, box->mid_z+extent);
                glEnd();

                glBegin(GL_LINE_LOOP);
                glVertex3f(box->mid_x+extent, box->mid_y+extent, box->mid_z+extent);
                glVertex3f(box->mid_x+extent, box->mid_y+extent, box->mid_z-extent);
                glVertex3f(box->mid_x+extent, box->mid_y-extent, box->mid_z-extent);
                glVertex3f(box->mid_x+extent, box->mid_y-extent, box->mid_z+extent);
                glVertex3f(box->mid_x+extent, box->mid_y+extent, box->mid_z+extent);
                glVertex3f(box->mid_x-extent, box->mid_y+extent, box->mid_z+extent);
                glVertex3f(box->mid_x-extent, box->mid_y+extent, box->mid_z-extent);
                glVertex3f(box->mid_x-extent, box->mid_y-extent, box->mid_z-extent);
                glVertex3f(box->mid_x-extent, box->mid_y-extent, box->mid_z+extent);
                glVertex3f(box->mid_x-extent, box->mid_y+extent, box->mid_z+extent);
                glEnd();
            }
            
            for(int i = 0 ; i < 8; ++i)
	    {
	      bool leaf = true;
	      if(box->child[i])
	      {
		leaf = false;
		num_boxes++;
		createBoxes(box->child[i],particles,extent*.5);
	      }
	      if(leaf)
	      {
		  int low = box->pArrayLow;
		  int high = box->pArrayHigh;
		  for (int p = low; p <= high; ++p)
		  {
		    GLfloat r = .005;
                    GLUquadric * qobj = gluNewQuadric();
                    glPushMatrix();
                    glTranslatef(particles[p].position[0], particles[p].position[1], particles[p].position[2]);
                    GLint slices = 8;
                    GLint stacks = 8;
                    gluSphere(qobj, r, slices, stacks);
                    glPopMatrix();
                    gluDeleteQuadric(qobj);
		    num_particles++;
		  }
	      }
		
	    }
            

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
            GLfloat x = 10.0 * GLfloat(width()) / height();
            glOrtho(-x, +x, -x, +x, -250.0, 200.0);

            glMatrixMode(GL_MODELVIEW);
            glPushMatrix();
            glLoadIdentity();
            glTranslatef(0.0, 0.0, -15.0);
            glScalef(m_scaling, m_scaling, m_scaling);

            glRotatef(m_rotationX, 1.0, 0.0, 0.0);
            glRotatef(m_rotationY, 0.0, 1.0, 0.0);
            glRotatef(m_rotationZ, 0.0, 0.0, 1.0);

            glEnable(GL_MULTISAMPLE);

            glCallList(m_glObject);

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

            QString label = "<h4 align=\"center\">Tree Info</h4>";
            label += "<p align=\"center\"><table width=\"100%\">";
            label += "<tr><td>Num Boxes: <td> "+ QString::number(num_boxes);
//             label += "<tr><td>Depth: <td> "+ QString::number(depth);
            label += "<tr><td>Num Particles: <td> "+ QString::number(num_particles);
//             label += "<tr><td>Particles per box: <td> "+ QString::number(max_particles);

            QTextDocument textDocument;
            textDocument.setDefaultStyleSheet("* { color: #FFEFEF }");
            textDocument.setHtml(label);
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
        int depth;
        int num_boxes;
        int num_particles;
        int max_particles;
};



#endif  // PARTICLE_RENDERER_HPP

