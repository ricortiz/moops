#ifndef MULTIPOLE_H
#define MULTIPOLE_H

#ifdef hpux
#define _HPUX_SOURCE 1
#endif
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#if defined(__STDC__) || defined(sgi) || defined(_AIX)
#include <unistd.h>
#else
#include <malloc.h>
#include <memory.h>
#endif
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
/*
  Define declarations for the N-body program.
*/
#define AbsoluteValue(x)  ((x) < 0 ? -(x) : (x))
#ifndef False
#define False  0
#endif
#define GravitationalConstant   6.673e-11
#define MaxPow  20
#define Max(x,y)  (((x) > (y)) ? (x) : (y))
#define MaxChildren  8
#define MaxNodesInteractiveField  189 
#define MaxNodesNearField  26  
#define MaxList3 500
#define MaxList4 2000
#define MaxList1 500
#define MaxColleagues 500
#define MaxList1Stack 500
#define MaxList2Stack 500
#define MaxTreeDepth  30
#define Epsilon 1.0e-5
#define MinRange  0.0
#define MaxRange  255.9999999
#define MEMORY_ALIGNMENT  4096
#define ALIGN_UP(x,size) ( ((size_t)x+(size-1))&(~(size-1)) )
#define Min(x,y)  (((x) < (y)) ? (x) : (y))
#ifndef RAND_MAX
#define RAND_MAX  32767
#endif
#define ArrayLookup(array,i,j,k)  \
  array[octree.ref_three[i]-octree.ref_two[(i)+(j)]+(k)]
#ifndef True
#define True  1
#endif

/*
  Structures.
*/
typedef struct _Particle
{
//   float position[3], force[3];
  float *position, *force;
} Particle;


typedef struct _Interval
{
    int low, high;

} Interval;

typedef struct _NodeInfo
{
    int low, high, level, subtree_size;
    unsigned int isParent;
    float mid_x, mid_y, mid_z;
    struct _NodeInfo *childInfo[MaxChildren];
    struct _Interval childInts[MaxChildren];
    
} NodeInfo;

typedef struct _Node
{
  struct _Node *parent, *child[MaxChildren], *Colleagues[MaxColleagues],
    *list1[MaxList1], *list2[MaxList3], *list3[MaxList3], *list4[MaxList4];
  unsigned long id,level;
  unsigned int DI_tag;
  float mid_x, mid_y, mid_z, *phi[4], *psi[4];
  int pArrayLow, pArrayHigh, *colleagueCount, *list1Count, *list2Count,
    *list4Count,*list3Count, isParent;
} Node;

typedef struct _Synch
{
    int flag;
    char pad[120]; 
} Synch;

typedef struct _GPU_Velocities
{
    float x,y,z;
} GPU_Velocities;

typedef struct _CPU_Velocities
{
    double x,y,z;
} CPU_Velocities;
typedef struct _Potential
{
    double x,y,z;
} Potential;

typedef struct _Field
{
    double field[3][4];
} Field;

typedef struct _Octree
{
    Node *root, *nodes, **DInodes, **leaf_DI_nodes;
    unsigned long depth;
    unsigned int *target_list, *number_IP, *interaction_pairs, 
        numParticles, number_DInodes, total_interaction_pairs, numLeafDInodes;
    int *counts;
    Node **near_neighbor;
    NodeInfo *rootInfo;
    FILE *output;
    float *GPU_Veloc;
    float *CPU_Veloc;
//     GPU_Velocities *GPU_Veloc;
//     CPU_Velocities *CPU_Veloc;
    Potential *potentials;
    Field *fields;    

    double DI_coeff, *scale;
    float edge_length[MaxTreeDepth], binomial[MaxPow+1][MaxPow+1],
        *ijk_binomial,*ijk_factorial,*PHIs[4], *PSIs[4];
    unsigned long ref_two[MaxPow+3], ref_three[MaxPow+3], precision,
        coefficients;
    Particle *bodies;

    int maxBodiesPerNode, numThreads, nodeInfoCoeff, nodeInfoDistribution
        ,*threadFailures, threads, treeRebuilds, DI_nodeTCount;

} Octree;
    /*
      Variable declarations.
    */
 
    Octree octree;
    /*
      Forward declarations.
    */
    unsigned int Binomial(int,int);
    void BuildColleagues(Node *);
    void BuildList1(Node *);
    void BuildList2(Node *);
    void BuildList3And4(Node *);
    void CheckColleagues(Node *, Node *);
    void CheckList3(Node *, Node *);
    void CheckList1(Node *, Node *);
    void UpSweep(Node *);
    void DownSweep(register Node *);
    void DirectViaExpansion(Node *, Node *);
    void PerformDirectInterations(Node *);
    void FindNeighbors(register Node *);
    void UpdateBodies(unsigned long,double);
    void UpdateBodiesWithGPU(unsigned long, double);
    void PrintBodies(int);
    int IsBorder(Node *, Node *);
    void MultipoleExpansion(Node *,Node *,register float *);
    void FMM(double);
    void AllPairs(unsigned long, double);
    void par_time(double *, double *, double *, double *);
    void PickVector(double rad, int ndim, float *x, float *y, float *z);
    void Plummer(int);
    double randdouble(void);
    void set_affinity(void);
    double wcTime(void);
    float Modulus(float a, float b);
    
    /*******************************************************************
     * 
     * For tree.c
     * 
     ******************************************************************/
    
    NodeInfo *AnalyzeTree(int, int);
    void constructTreeInfo(int, int);
    void CreateSubtree(Node *, NodeInfo *, int);
    void CreateOctree(unsigned long, int, double, double);
    Node *InitializeNode(int, int ,Node *,double,double,double, int, int, int, int);
    void ParticleSwap(int, int);
    int PartitionBodies(int,int,float, int);
    void PartitionInfo(NodeInfo *);    
    void subdivide(int, int, float, float, float, Interval *);
    void threadSetup();
    void treeVarSetup(int,int, double, double);
    void CountDINodes(Node *);
    void constructDINodes(unsigned long, double);
    void ConstructDINodesRec(Node *);
    void ParPrefixSum(unsigned int *, unsigned int);
    int ReSort(Node *node, NodeInfo *info);
    void PushDownDirectWork(Node*);
    void FreeOctree();
    void RebuildTree(unsigned long number_particles, int precision, double maximum_extent, double minimum_extent);

    /*******************************************************************
     *
     * For expansions.c
     *
     ******************************************************************/

     void ApplyLocalExpansion(Node *);
     void CalculateDirectInteractions(Node *, Node **, int);
     void FormOuterExpansion(Node *);
     void FormLocalExpansion(Node *,Node *, float *);
     void ShiftFromChildToParent(Node *);
     void ShiftFromParentToLeaf(Node *);
     void PerformDirectWithinNode(Node *);
     void OuterToInner(Node *source, Node *target, float *target_psi[4], int isLeaf);
     void DownShift(Node *, int);
     void BuildGPUArrays();
     void ComputeVelocityDirect(int target, int source);
     

     /******************************************************************
      *
      * From libkernel.a
      *
      *****************************************************************/     
    void gpuVelocitiesEval( const unsigned int NUM_BODIES, 
                        const unsigned int NUM_LEAF_NODES, 
                        const unsigned int TOTAL_NUM_SOURCES,
                        float * positions, 
                        float * gpuVelocities, 
                        unsigned int * target_list, 
                        unsigned int * num_interaction_pairs, 
                        unsigned int * interaction_pairs, float delta);

    void gpuGetVelocities();


#endif
