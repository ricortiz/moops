
#ifndef TYPE_BINDER_HPP
#define TYPE_BINDER_HPP

#ifdef __INTEL_COMPILER
#pragma warning(disable:193 383 444 981 1572 2259)
#endif

#include <algorithm>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <omp.h>
#include <queue>
#include <stack>
#include <string>
#include <utility>
#include <vector>
#include "vec.h"                                                //!< My vector type with operator overloading
#if PAPI
#include <papi.h>
#endif
#if QUARK
#include "quark.h"
#endif

template<typename real_type>
struct default_math_type
{
    typedef unsigned    int_type;
    typedef real_type   real_type;
    typedef std::complex<real_type> complex_type;
    typedef vec<3, real_type> vector3_type;
    typedef vec<4, real_type> vector4_type;
    typedef vec<8, real_type> vector8_type;
};

//! Structure for pthread based trace
struct Trace
{
    pthread_t thread;
    double    begin;
    double    end;
    int       color;
};

enum Equation                                                   //!< Equation type enumeration
{
    Laplace,                                                      //!< Laplace potential + force
    VanDerWaals,                                                  //!< Van der Walls potential + force
    Stokes                                                        //!< Stokes flow
};

//! Structure of source bodies (stuff to send)
template<typename types>
struct SourceBody
{
    typedef typename types::int_type int_type;
    typedef typename types::vector3_type vector_3_type;
    typedef typename types::real_type real_type;
    int         IBODY;                                            //!< Initial body numbering for sorting back
    int         IPROC;                                            //!< Initial process numbering for sending back
    int_type      ICELL;                                            //!< Cell index
    vector_3_type        X;                                                //!< Position
    real_type        SRC;                                              //!< Scalar source values
    SourceBody() : IBODY(0), IPROC(0), ICELL(0), X(0), SRC(0) {}       //!< Constructor
};

//! Structure of bodies
template<typename types>
struct TargetBody : public SourceBody<types>
{
    vec<4, real> TRG;                                             //!< Scalar+vector target values
    bool operator<(const Body &rhs) const                         //!< Overload operator for comparing body index
    {
        return this->IBODY < rhs.IBODY;                             //!< Comparison function for body index
    }                                                             //!< End operator overload
    Body() : TRG(0) {}                                            //!< Constructor
};

//! Linked list of leafs (only used in fast/topdown.h)
template<typename types>
struct Leaf
{
    typedef typename types::vector3_type vector_3_type;
    int I;                                                        //!< Unique index for every leaf
    vector_3_type X;                                                       //!< Coordinate of leaf
    Leaf *NEXT;                                                   //!< Pointer to next leaf
    Leaf() : I(0), X(0), NEXT(NULL) {}                            //!< Constructor
    ~Leaf() {}                                                    //!< Destructor
//! Copy constructor
    Leaf(const Leaf &leaf) : I(0), X(0), NEXT(NULL)
    {
        I    = leaf.I;                                              // Copy I
        X    = leaf.X;                                              // Copy X
        NEXT = leaf.NEXT;                                           // Copy NEXT
    }
//! Overload assignment
    Leaf &operator=(const Leaf &leaf)
    {
        I    = leaf.I;                                              // Copy I
        X    = leaf.X;                                              // Copy X
        NEXT = leaf.NEXT;                                           // Copy NEXT
        return *this;
    }
};

//! Structure of nodes (only used in fast/topdown.h)
template<typename types>
struct Node
{
    typedef typename types::vector3_type vector_3_type;
    typedef typename types::leaf_type leaf_type;
    bool NOCHILD;                                                 //!< Flag for twig nodes
    int  LEVEL;                                                   //!< Level in the tree structure
    int  NLEAF;                                                   //!< Number of descendant leafs
    vec<8, int> CHILD;                                            //!< Index of child node
    vector_3_type X;                                                       //!< Coordinate at center
    leaf_type *LEAF;                                                   //!< Pointer to first leaf
    Node() : NOCHILD(true), LEVEL(0), NLEAF(0), CHILD(-1), X(0), LEAF(NULL) {}//!< Constructor
    ~Node() {}                                                    //!< Destructor
//! Copy constructor
    Node(const Node &node) : NOCHILD(true), LEVEL(0), NLEAF(0), CHILD(-1), X(0), LEAF(NULL)
    {
        NOCHILD = node.NOCHILD;                                     // Copy NOCHILD
        LEVEL   = node.LEVEL;                                       // Copy LEVEL
        NLEAF   = node.NLEAF;                                       // Copy NLEAF
        CHILD   = node.CHILD;                                       // Copy CHILD
        X       = node.X;                                           // Copy X
        LEAF    = node.LEAF;                                        // Copy LEAF
    }
//! Overload assignment
    Node &operator=(const Node &node)
    {
        NOCHILD = node.NOCHILD;                                     // Copy NOCHILD
        LEVEL   = node.LEVEL;                                       // Copy LEVEL
        NLEAF   = node.NLEAF;                                       // Copy NLEAF
        CHILD   = node.CHILD;                                       // Copy CHILD
        X       = node.X;                                           // Copy X
        LEAF    = node.LEAF;                                        // Copy LEAF
        return *this;
    }
};

//! Structure of source cells (stuff to send)
template<typename types>
struct SourceCell
{
    typedef typename types::int_type int_type;
    typedef typename types::Mset Mset;
    int_type ICELL;                                                 //!< Cell index
    Mset   M;                                                     //!< Multipole coefficients
};

//! Structure of cells
template<typename types>
struct Cell
{
    typedef typename types::int_type int_type;
    typedef typename types::real_type real_type;
    typedef typename types::vector3_type vector3_type;
    typedef typename types::B_iter body_iterator;
    typedef typename types::Mset multipole_expansion;
    typedef typename types::Lset local_expansion;
    int_type ICELL;                                                 //!< Cell index
    int    NCHILD;                                                //!< Number of child cells
    int    NCLEAF;                                                //!< Number of child leafs
    int    NDLEAF;                                                //!< Number of descendant leafs
    int    PARENT;                                                //!< Iterator offset of parent cell
    int    CHILD;                                                 //!< Iterator offset of child cells
    int    ILEAF;                                                 //!< Iterator offset of first leaf
    body_iterator LEAF;                                                  //!< Iterator of first leaf
    vector3_type   X;                                                     //!< Cell center
    real_type   R;                                                     //!< Cell radius
    real_type   RMAX;                                                  //!< Max cell radius
    real_type   RCRIT;                                                 //!< Critical cell radius
    multipole_expansion  M;                                                     //!< Multipole coefficients
    local_expansion   L;                                                     //!< Local coefficients
    Cell() : ICELL(0), NCHILD(0), NCLEAF(0), NDLEAF(0), PARENT(0), CHILD(0),
            LEAF(), X(0), R(0), RMAX(0), RCRIT(0), M(0), L(0) {} //!< Constructor
};

//! Structure for Ewald summation
template<typename types>
struct Ewald
{
    typedef typename types::real_type real_type;
    typedef typename types::vector3_type vector3_type;
    vector3_type K;                                                       //!< 3-D wave number vector
    real_type REAL;                                                    //!< Real part of wave
    real_type IMAG;                                                    //!< Imaginary part of wave
    Ewald() : K(0), REAL(0), IMAG(0) {}                           //!< Constructor
};

template <typename math_types, int _P>
struct TypeBinder
{
    typedef TypeBinder<math_type, _P> types;
    enum
    {
        P = _P,                                         //!< Order of expansions
        MTERM = P * (P + 1) * (P + 2) / 6,              //!< Number of Cartesian mutlipole terms
        LTERM = (P + 1) * (P + 2) * (P + 3) / 6,        //!< Number of Cartesian local terms
        NTERM = P * (P + 1) / 2,                        //!< Number of Spherical multipole/local terms
        NCRIT = 10,                                     //!< Number of bodies per cell
        MAXBODY = 2000000,                              //!< Maximum number of bodies per GPU kernel
        MAXCELL = 10000000,                             //!< Maximum number of bodies/coefs in cell per GPU kernel
        GPUS = 2,                                       //!< Number of GPUs per node
        THREADS  = 64,                                  //!< Number of threads per thread-block
        PTHREADS = 4                                    //!< Number of pthreads in quark
    };
    typedef typename math_types::int_type       int_type;
    typedef typename math_types::real_type      real_type;      //!< Real number type
    typedef typename math_types::complex_type   complex_type;   //!< Complex number type
    typedef typename math_types::vector3_type   vector3_type;   //!< 3-D vector type
    typedef std::vector<int_type>               int_array_type;    //!< Vector of big integer types
    typedef std::map<pthread_t, double>         ThreadTrace;    //!< Map of pthread id to traced value
    typedef std::map<pthread_t, int>            ThreadMap;      //!< Map of pthread id to thread id
    typedef std::queue<Trace>                   Traces;          //!< Queue of traces
    typedef std::map<std::string, double>       Timer;          //!< Map of timer event name to timed value
    typedef std::map<std::string, double>::iterator TI_iter;        //!< Iterator for timer event name map

    typedef SourceBody<types>           JBody;
    typedef std::vector<JBody>             JBodies;                 //!< Vector of source bodies
    typedef std::vector<JBody>::iterator   JB_iter;                 //!< Iterator for source body vector

    typedef TargetBody<types>           Body;
    typedef std::vector<Body>              Bodies;                  //!< Vector of bodies
    typedef std::vector<Body>::iterator    B_iter;                  //!< Iterator for body vector

    typedef Leaf<types>                 leaf_type;
    typedef std::vector<leaf_type>              Leafs;                   //!< Vector of leafs
    typedef std::vector<leaf_type>::iterator    L_iter;                  //!< Iterator for leaf vector

    typedef Node<types>                 node_type;
    typedef std::vector<node_type>              Nodes;                   //!< Vector of nodes
    typedef std::vector<node_type>::iterator    N_iter;                  //!< Iterator for node vector

    typedef SourceCell<types>   JCell;
    typedef std::vector<JCell>             JCells;                  //!< Vector of source cells
    typedef std::vector<JCell>::iterator   JC_iter;                 //!< Iterator for source cell vector

    typedef Cell<types>         cell_type;
    typedef std::vector<cell_type>              Cells;                   //!< Vector of cells
    typedef std::vector<cell_type>::iterator    C_iter;                  //!< Iterator for cell vector

    typedef Ewald<types>        ewald_type;
    typedef std::vector<ewald_type>             Ewalds;                  //!< Vector of Ewald summation types
    typedef std::vector<ewald_type>::iterator   E_iter;                  //!< Iterator for Ewald summation types

    typedef std::queue<C_iter>             CellQueue;               //!< Queue of cell iterators
    typedef std::stack<C_iter>             CellStack;               //!< Stack of cell iterators
    typedef std::pair<C_iter, C_iter>       Pair;                   //!< Pair of interacting cells
    typedef std::queue<Pair>               PairQueue;               //!< Queue of interacting cell pairs
    typedef std::stack<Pair>               PairStack;               //!< Stack of interacting cell pairs
    typedef std::list<C_iter>              List;                    //!< Interaction list
    typedef std::list<C_iter>::iterator    LC_iter;                 //!< Iterator for interaction list
    typedef std::vector<List>              Lists;                   //!< Vector of interaction lists
    typedef std::map<C_iter, int>           Map;                    //!< Map of interaction lists
    typedef std::map<C_iter, int>::iterator MC_iter;                //!< Iterator for interation list map
    typedef std::vector<Map>               Maps;                    //!< Vector of map of interaction lists

#if Cartesian
    typedef vec<MTERM, real_type>                        Mset;     //!< Multipole coefficient type for Cartesian
    typedef vec<LTERM, real_type>                        Lset;     //!< Local coefficient type for Cartesian
#elif Spherical
    typedef vec<NTERM, complex_type>                     Mset;     //!< Multipole coefficient type for spherical
    typedef vec<NTERM, complex_type>                     Lset;     //!< Local coefficient type for spherical
#endif

#ifndef KERNEL
    int MPIRANK    = 0;                                             //!< MPI comm rank
    int MPISIZE    = 1;                                             //!< MPI comm size
    int DEVICE     = 0;                                             //!< GPU device ID
    int IMAGES     = 0;                                             //!< Number of periodic image sublevels
    real_type THETA     = .5;                                            //!< Multipole acceptance criteria
    vector3_type Xperiodic = 0;                                             //!< Coordinate offset of periodic image
#if PAPI
    int PAPIEVENT  = PAPI_NULL;                                     //!< PAPI event handle
#endif
#else
    extern int MPIRANK;                                             //!< MPI comm rank
    extern int MPISIZE;                                             //!< MPI comm size
    extern int DEVICE;                                              //!< GPU device ID
    extern int IMAGES;                                              //!< Number of periodic image sublevels
    extern real_type THETA;                                         //!< Multipole acceptance criteria
    extern vector3_type Xperiodic;                                  //!< Coordinate offset of periodic image
#if PAPI
    extern int PAPIEVENT;                                           //!< PAPI event handle
#endif
#endif
    const real_type CLET     = 2;                                        //!< LET opening critetia
    const real_type EPS      = 1e-6;                                     //!< Single precision epsilon
    const real_type EPS2     = 0;                                        //!< Softening parameter (squared)
    const real_type R2MIN    = 0.25;                                     //!< Minimum value for L-J R^2
    const real_type R2MAX    = 64;                                       //!< Maximum value for L-J R^2

};

#endif
