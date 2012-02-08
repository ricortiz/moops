/*
  Include declarations.
*/
#define _GNU_SOURCE
#include <omp.h>
#include "hybridfmm.h"
#include <sys/time.h>
#ifndef __INTEL_COMPILER
#include <math.h>
#else
#include "patched_math.h"
#endif
#include <assert.h>
#include <stdio.h>
#include <sched.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <linux/unistd.h>

//@TODO: easy way than using a stack???? re-do like the in build list3and4

/*
  Global variables.
*/
FILE * gpuOut;
double upTime, downTime;
#pragma omp threadprivate(upTime, downTime)
/****************************************************************
 *
 * returns binomial coefficient
 *
 ****************************************************************/
unsigned int Binomial(int n, int k)
{
    int i, j;
    unsigned int binomial;

    if (n < k)
    {
        return 0;
    }
    binomial = 1;
    k = Min(k, n - k);
    j = n;
    for (i = 1; i <= k; i++)
    {
        binomial = (unsigned long)((binomial * (float) j / i) + 0.5);
        j--;
    }
    return binomial;
}

/***********************************************************************
 *
 * UpSweep( Node *node)
 * @brief - Carries out the upsweep phase of FMM. Upsweet is achieved
 * by first recursing to bottom and then returning upwards
 *
 **********************************************************************/
void UpSweep(Node *node)
{
    int id;
    double stime;

    if (node->isParent)
    {
        for (id = 0; id < MaxChildren; id++)
        {
            if (node->child[id]->pArrayLow != -1)
            {
#pragma omp task firstprivate(id)
                UpSweep(node->child[id]);
            }
        }
#pragma omp taskwait
    }
    stime = wcTime();
    if (!node->isParent)
    {
        FormOuterExpansion(node);
    }
    if (node->isParent && (node->level > 0))
    {
        ShiftFromChildToParent(node);
    }
    upTime += wcTime() - stime;
}

/***********************************************************************
 *
 * DownSweep(Node *node)
 *
 **********************************************************************/
void DownSweep(Node *node)
{
    int id, i;
    double stime;

    stime = wcTime();
    if (!node->isParent)
    {
        //shift in list 2
        for (i = 0; i < *(node->list2Count); i++)
        {
            OuterToInner(node->list2[i], node, node->psi, True);
        }
        DownShift(node, True);
        for (i = 0;i < *node->list3Count;i++)
        {
            if (node->list3[i]->isParent)
                OuterToInner(node->list3[i], node, node->psi, True);
        }
        ApplyLocalExpansion(node);
    }
    else if (node->level > 1)
    {
        for (i = 0; i < *(node->list2Count); i++)
        {
            OuterToInner(node->list2[i], node, node->psi, False);
        }
        DownShift(node, False);
    }
    downTime += wcTime() - stime;
    //recurse
    if (node->isParent)
    {
        for (id = 0; id < MaxChildren; id++)
        {
            if (node->child[id]->pArrayLow != -1)
            {
#pragma omp task firstprivate(id)
                DownSweep(node->child[id]);
            }
        }
#pragma omp taskwait
    }

}

void CheckColleagues(Node *node, Node *candidate)
{
    int j;
    if ((candidate->level < node->level) && candidate->isParent)
    {
        for (j = 0; j < MaxChildren; j++)
        {
            if ((candidate->child[j] != NULL) && (IsBorder(node, candidate->child[j]) == 1))
            {
                CheckColleagues(node, candidate->child[j]);
            }
        }
    }
    else
    {
        if (IsBorder(node, candidate) == 1)
        {
            node->Colleagues[*node->colleagueCount] = candidate;
            *node->colleagueCount = *node->colleagueCount + 1;
        }
    }
}
/**
void BuildColleagues(Node *node){
    int i,j,k,count;
    Node **colleaguesList, **srcColleaguesList;
    if(node == octree.root)
    {
        for(i = 0; i < MaxChildren; i++)
        {
            count = 0;
            if(octree.root->child[i]!=NULL)
            {
                colleaguesList = octree.root->child[i]->Colleagues;
                for(j=0; j<MaxChildren; j++)
                {
                    if((octree.root->child[j]!=NULL)
                    &&(octree.root->child[j]!=octree.root->child[i])){
                        colleaguesList[count++] = octree.root->child[j];
                    }
                }
                for(k=0;k<MaxChildren;k++){
                    if(octree.root->child[i]->child[k]!=NULL){
                        BuildColleagues(octree.root->child[i]->child[k]);
                    }
                }
            }
        }
    }
    else
    {
        for(i=0;i<MaxColleagues;i++){
            if(node->parent->Colleagues[i]==NULL){
                break;
            }
            else if((node->parent->Colleagues[i]->isParent)
                    && (IsBorder(node, node->parent->Colleagues[i])==1)
                    && (node->parent->Colleagues[i]->level < node->level))
            {
                for(j=0;j<MaxChildren;j++){
                    if((node->parent->Colleagues[i]->child[j]!=NULL)
                        &&(IsBorder(node, node->parent->Colleagues[i]->child[j])==1))
                    {
                        node->Colleagues[*node->colleagueCount] = node->parent->Colleagues[i]->child[j];
                        *node->colleagueCount = *node->colleagueCount + 1;
                    }
                }
            }
            else if(IsBorder(node, node->parent->Colleagues[i])==1){
                 node->Colleagues[*node->colleagueCount] = node->parent->Colleagues[i];
                *node->colleagueCount = *node->colleagueCount + 1;
            }
        }
        for(i=0;i<MaxChildren;i++){
            if((node->parent->child[i]!=NULL)&&(node->parent->child[i]!=node))
            {
                node->Colleagues[*node->colleagueCount] = node->parent->child[i];
                *node->colleagueCount = *node->colleagueCount + 1;
            }
        }
    }
    if((node->level<(octree.depth-1))&&(node->isParent)){
        for(k=0;k<MaxChildren;k++){
            if(node->child[k]!=NULL){
                //#pragma omp task firstprivate(k)
                BuildColleagues(node->child[k]);
            }
        }
        //#pragma omp taskwait
    }
}*/

/***********************************************************************
 *
 * Builds the colleagues list for each node. This list is the list
 * of nodes that border a given node. This includes those nodes which
 * are  adjacent to the given node as well as those that only shared a
 * point. Nodes in this list are not necessarily on the same level as
 * the given node nor are they necessarily childless. This list is used
 * to build Lists 1 through 4.
 *
 **********************************************************************/
void BuildColleagues(Node * node)
{
    int i, j, k, count;
    Node **colleaguesList, **srcColleaguesList;

    if (node == octree.root)
    {
        for (i = 0; i < MaxChildren; i++)
        {
            count = 0;
            if (octree.root->child[i] != NULL)
            {
                colleaguesList = octree.root->child[i]->Colleagues;
                for (j = 0; j < MaxChildren; j++)
                {
                    if ((octree.root->child[j] != NULL)
                            && (octree.root->child[j] != octree.root->child[i]))
                    {
                        colleaguesList[count++] = octree.root->child[j];
                    }
                }
                for (k = 0;k < MaxChildren;k++)
                {
                    if (octree.root->child[i]->child[k] != NULL)
                    {
                        BuildColleagues(octree.root->child[i]->child[k]);
                    }
                }
            }
        }
    }
    else
    {
        Node *stack[MaxList1Stack], *temp;
        int top;
        top = -1;
        for (i = 0; i < (MaxList1Stack); i++)
        {
            stack[i] = NULL;
        }

        srcColleaguesList = node->parent->Colleagues;
        colleaguesList = node->Colleagues;
        count = 0;
        for (i = 0; i < (MaxList3); i++) //check parents colleagues
        {
            if (srcColleaguesList[i] != NULL)
            {
                //put any node that borders @node on the stack
                //before adding to node->BorderList we will check to
                //see if it can be subdivided more
                if (IsBorder(node, srcColleaguesList[i]) == 1)
                {
                    stack[++top] = srcColleaguesList[i];
                }
            }
        }
        while (top != -1)
        {
            //if it has children and is at earlier level than node
            if (stack[top]->isParent && (stack[top]->level < node->level))
            {
                //take off stack since only interested in its children now
                temp = stack[top--];
                for (j = 0;j < MaxChildren;j++)
                {
                    if (temp->child[j] != NULL)
                    {
                        if (IsBorder(node, temp->child[j]) == 1)//only add if borders
                        {
                            stack[++top] = temp->child[j];
                        }
                    }
                }
            }
            else  //put in colleagues list
            {
                //assert(stack[top]->level<=node->level);
                colleaguesList[count++] = stack[top];
                stack[top--] = NULL;
            }

        }
        for (i = 0; i < MaxChildren; i++)
        {
            if ((node->parent->child[i] != NULL) && (node->parent->child[i] != node))
            {
                colleaguesList[count++] = node->parent->child[i];
            }
        }
        if ((node->level < (octree.depth - 1)) && (node->isParent))
        {
            for (k = 0;k < MaxChildren;k++)
            {
#pragma omp task firstprivate(k)
                BuildColleagues(node->child[k]);
            }
#pragma omp taskwait
        }
    }
}

/***********************************************************************
 *
 * Checks if candidate is a parent, if so for each child of candidate
 * adjacent to node, call CheckList1 on it. If not parent, then must be
 * adjacent to node so add to list
 *
 **********************************************************************/
void CheckList1(Node *node, Node *candidate)
{
    int j;
    if (candidate->isParent)
    {
        for (j = 0; j < MaxChildren; j++)
        {
            if (IsBorder(node, candidate->child[j]) == 1)
            {
                CheckList1(node, candidate->child[j]);
            }
        }
    }
    else
    {
        node->list1[*node->list1Count] = candidate;
        *node->list1Count = *node->list1Count + 1;
    }
}

/*********************************************************************
 *
 * Builds List1 for a node and its children
 *
 ********************************************************************/
void BuildList1(Node *node)
{
    int i;
    if (!node->isParent)
    {
        //check each colleague
        for (i = 0; i < MaxColleagues; i++)
        {
            //then seen all colleagues
            if (node->Colleagues[i] == NULL)
            {
                break;
            }
            else
            {
                CheckList1(node, node->Colleagues[i]);
            }
        }
    }

    if (node->isParent)
    {
        for (i = 0; i < MaxChildren; i++)
        {
#pragma omp task firstprivate(i)
            BuildList1(node->child[i]);
        }
#pragma omp taskwait
    }
}
/**********************************************************************
 *
 * Builds List two for a node and the recursively descends
 *
 *********************************************************************/
void BuildList2(Node *node)
{
    int i, j;

    for (i = 0;i < MaxColleagues;i++)
    {
        if (node->parent->Colleagues[i] == NULL)
        {
            break;
        }
        else if (node->parent->level == node->parent->Colleagues[i]->level) //then true colleague
        {
            if (node->parent->Colleagues[i]->isParent)
            {
                for (j = 0;j < MaxChildren;j++)
                {
                    if ((node->parent->Colleagues[i]->child[j] != NULL)
                            && (IsBorder(node, node->parent->Colleagues[i]->child[j]) != 1))
                    {
                        node->list2[*node->list2Count] = node->parent->Colleagues[i]->child[j];
                        *node->list2Count = *node->list2Count + 1;
                    }
                }
            }
        }
    }
    //descend
    if (node->isParent)
    {
        for (i = 0; i < MaxChildren; i++)
        {
#pragma omp task firstprivate(i)
            BuildList2(node->child[i]);
        }
#pragma omp taskwait
    }
}

void CheckList3(Node * node, Node * candidate)
{
    int j;
    if (IsBorder(node, candidate) != 1) //if not adjacent
    {
        //put in list3 and by result list4
        node->list3[*node->list3Count] = candidate;
        *node->list3Count = *node->list3Count + 1;
        candidate->list4[*candidate->list4Count] = node;
        *candidate->list4Count = *candidate->list4Count + 1;
    }
    else if (candidate->isParent)
    {
        for (j = 0; j < MaxChildren; j++)
        {
            CheckList3(node, candidate->child[j]);
        }
    }

}

/***********************************************************************
 *
 * BuildList3And4 builds 3&4
 *
 **********************************************************************/
void BuildList3And4(Node * node)
{
    int i, j;

    if (!node->isParent)
    {
        //for each colleague
        for (i = 0; i < MaxColleagues; i++)
        {
            //then seen all colleagues
            if (node->Colleagues[i] == NULL)
            {
                break;
            }
            else if (node->Colleagues[i]->isParent) //each child of colleague is potential list3 node
            {
                for (j = 0; j < MaxChildren; j++)
                {
                    CheckList3(node, node->Colleagues[i]->child[j]);
                }
            }
        }
    }

    if (node->isParent)//then has at least one child
    {
        for (i = 0; i < MaxChildren; i++)
        {
            //#pragma omp task firstprivate(i)
            BuildList3And4(node->child[i]);
        }
        //#pragma omp taskwait
    }
}

/***********************************************************************
 *
 * IsBorder(Node * constructListFor, Node * toCheck) checks if toCheck
 * borders constructListFor and returns 1 if so, 0 otherwise.
 *
 **********************************************************************/
int IsBorder(Node * constructListFor, Node * toCheck)
{

    float toCheckMidX, toCheckMidY, toCheckMidZ, constructMidX, constructMidY, constructMidZ,
    minusCheckX, plusCheckX, minusCheckY, plusCheckY, minusCheckZ, plusCheckZ, minusConstX,
    plusConstX, minusConstY, plusConstY, minusConstZ, plusConstZ;
    double toCheckBisect, constructBisect;

    int toReturn;

    toCheckMidX = toCheck->mid_x;
    toCheckMidY = toCheck->mid_y;
    toCheckMidZ = toCheck->mid_z;
    toCheckBisect = octree.edge_length[toCheck->level] * 0.5;

    minusCheckX = toCheckMidX - toCheckBisect;
    plusCheckX = toCheckMidX + toCheckBisect;
    minusCheckY = toCheckMidY - toCheckBisect;
    plusCheckY = toCheckMidY + toCheckBisect;
    minusCheckZ = toCheckMidZ - toCheckBisect;
    plusCheckZ = toCheckMidZ + toCheckBisect;

    constructMidX = constructListFor->mid_x;
    constructMidY = constructListFor->mid_y;
    constructMidZ = constructListFor->mid_z;
    constructBisect = octree.edge_length[constructListFor->level] * 0.5;

    minusConstX = constructMidX - constructBisect;
    plusConstX = constructMidX + constructBisect;
    minusConstY = constructMidY - constructBisect;
    plusConstY = constructMidY + constructBisect;
    minusConstZ = constructMidZ - constructBisect;
    plusConstZ = constructMidZ + constructBisect;

    if ((minusConstX == plusCheckX) || (plusConstX == minusCheckX))
    {
        if ((plusConstY >= minusCheckY) && (minusConstY <= plusCheckY)
                && (minusConstZ <= plusCheckZ) && (plusConstZ >= minusCheckZ))
        {
            return 1;
        }
    }
    else if ((minusConstZ == plusCheckZ) || (plusConstZ == minusCheckZ))
    {
        if ((plusConstX >= minusCheckX) && (minusConstX <= plusCheckX)
                && (minusConstY <= plusCheckY) && (plusConstY >= minusCheckY))
        {
            return 1;
        }
    }
    else if ((minusConstY == plusCheckY) || (plusConstY == minusCheckY))
    {
        if ((plusConstX >= minusCheckX) && (minusConstX <= plusCheckX)
                && (plusConstZ >= minusCheckZ) && (minusConstZ <= plusCheckZ))
        {
            return 1;
        }
    }
    return 0;
}

/***********************************************************************
 *
 * Does pairwise interactions for each particle in node with each
 * particle in each box in node's @list
 *
 **********************************************************************/
void CalculateDirectInteractions(Node *node, Node **list, int max, double delta)
{

    Node *list_node;
    int i, low, high, pNum, qNum, list_Low, list_High;

    low = node->pArrayLow;
    high = node->pArrayHigh;


    //for every node in list
    for (i = 0; i < max; i++)
    {
        list_node = list[i];
        {
            list_Low = list_node->pArrayLow;
            list_High = list_node->pArrayHigh;
            {
                if ((!list_node->isParent) && (list_Low != -1))
                {
                    //for every particle in this node
                    for (pNum = low; pNum <= high; pNum++)
                    {
                        //for every particle in list node
                        for (qNum = list_Low; qNum <= list_High; qNum++)
                        {
                            ComputeVelocityDirect(pNum, qNum, delta);
                        }
                    }
                }
            }
        }
    }
}

/***********************************************************************
 *
 * DirectViaExpansion(Node *target, Node *source)
 * @brief - supposed to replace direct interaction between the particles
 * of two nodes with a particle - expansion interaction but not working
 * yet so just set to allways do all pairs method
 *
 **********************************************************************/
void DirectViaExpansion(Node *target, Node *source, double delta)
{
    int i, j;

    for (i = target->pArrayLow; i <= target->pArrayHigh; i++)
    {
        for (j = source->pArrayLow; j <= source->pArrayHigh; j++)
        {
            ComputeVelocityDirect(i, j, delta);
        }
    }
    /*
    double dx, dy, dz, sTime;
    Particle *p;
    int i, j, k,low, high, pNum;
    float x_power[MaxPow+1], y_power[MaxPow+1], z_power[MaxPow+1], *psi_0,
        *psi_1, *psi_2, *psi_3,dx_power[MaxPow+1],dy_power[MaxPow+1],
        dz_power[MaxPow+1];

    //sTime = wcTime();
    float *temp_psi_0, *temp_psi_1, *temp_psi_2, *temp_psi_3;

    temp_psi_0 = (float *) malloc (octree.coefficients*sizeof(float));
    temp_psi_1 = (float *) malloc (octree.coefficients*sizeof(float));
    temp_psi_2 = (float *) malloc (octree.coefficients*sizeof(float));
    temp_psi_3 = (float *) malloc (octree.coefficients*sizeof(float));
    float *temp_psi[4] = {temp_psi_0,temp_psi_1,temp_psi_2,temp_psi_3};
    float tempPotentials[3] = {0};
    float tempField[3][4] = {{0}};

    memset(temp_psi_0, 0.0, octree.coefficients*sizeof(float));
    memset(temp_psi_1, 0.0, octree.coefficients*sizeof(float));
    memset(temp_psi_2, 0.0, octree.coefficients*sizeof(float));
    memset(temp_psi_3, 0.0, octree.coefficients*sizeof(float));


    low = target->pArrayLow;
    high = target->pArrayHigh;
    OuterToInner(source, target, temp_psi,True);

    for (pNum = low; pNum<= high; pNum++)
    {
        p = &octree.bodies[pNum];
        dx=p->position[0]-target->mid_x + Epsilon;
        dy=p->position[1]-target->mid_y + Epsilon;
        dz=p->position[2]-target->mid_z + Epsilon;
        for (i=0; i <= octree.precision; i++)
        {
            x_power[i]=pow(dx,(double) i);
            y_power[i]=pow(dy,(double) i);
            z_power[i]=pow(dz,(double) i);

            dx_power[i] = i*pow(dx,(double) i-1.0);
            dy_power[i] = i*pow(dy,(double) i-1.0);
            dz_power[i] = i*pow(dz,(double) i-1.0);
        }
        psi_0 = temp_psi_0;
        psi_1 = temp_psi_1;
        psi_2 = temp_psi_2;
        psi_3 = temp_psi_3;
        //No need for temps if it is safe to use existing, need to figure out if ok
        for (i=0; i <= octree.precision; i++){
            for (j=0; j <= (octree.precision-i); j++){
                for (k=0; k <= (octree.precision-i-j); k++)
                {
                    tempPotentials[0]+=(*psi_0)*x_power[i]*y_power[j]*z_power[k];
                    tempPotentials[1]+=(*psi_1)*x_power[i]*y_power[j]*z_power[k];
                    tempPotentials[2]+=(*psi_2)*x_power[i]*y_power[j]*z_power[k];

                    tempField[0][0] += (*psi_0)*dx_power[i]*y_power[j]*z_power[k];
                    tempField[1][0] += (*psi_0)*x_power[i]*dy_power[j]*z_power[k];
                    tempField[2][0] += (*psi_0)*x_power[i]*y_power[j]*dz_power[k];

                    tempField[0][1] += (*psi_1)*dx_power[i]*y_power[j]*z_power[k];
                    tempField[1][1] += (*psi_1)*x_power[i]*dy_power[j]*z_power[k];
                    tempField[2][1] += (*psi_1)*x_power[i]*y_power[j]*dz_power[k];

                    tempField[0][2] += (*psi_2)*dx_power[i]*y_power[j]*z_power[k];
                    tempField[1][2] += (*psi_2)*x_power[i]*dy_power[j]*z_power[k];
                    tempField[2][2] += (*psi_2)*x_power[i]*y_power[j]*dz_power[k];

                    tempField[0][3] += (*psi_3)*dx_power[i]*y_power[j]*z_power[k];
                    tempField[1][3] += (*psi_3)*x_power[i]*dy_power[j]*z_power[k];
                    tempField[2][3] += (*psi_3)*x_power[i]*y_power[j]*dz_power[k];

                    psi_0++;
                    psi_1++;
                    psi_2++;
                    psi_3++;
                }
            }
        }
        octree.CPU_Veloc[pNum].x += tempPotentials[0] - p->position[0]*tempField[0][0] - p->position[1]*tempField[0][1] - p->position[2]*tempField[0][2] + tempField[0][3];
        octree.CPU_Veloc[pNum].y += tempPotentials[1] - p->position[0]*tempField[1][0] - p->position[1]*tempField[1][1] - p->position[2]*tempField[1][2] + tempField[1][3];
        octree.CPU_Veloc[pNum].z += tempPotentials[2] - p->position[0]*tempField[2][0] - p->position[1]*tempField[2][1] - p->position[2]*tempField[2][2] + tempField[2][3];
    }
    //fprintf(octree.output, "%f, %i,%i\n", wcTime()- sTime, high - low + 1, source->pArrayHigh-source->pArrayLow+1);
    * */
}

/***********************************************************************
 *
 * FMM
 * @brief - carries out one FMM iteration.
 *
 **********************************************************************/
void FMM(double dt)
{
    double dTime, psiTime, phiTime;
    int i;

#pragma omp parallel
    {
        upTime = downTime = 0.0;
#pragma omp single
        {

            dTime = wcTime();
            gpuVelocitiesEval(octree.numParticles, octree.numLeafDInodes,
                              octree.total_interaction_pairs, (float *) octree.bodies,
                              (float *) octree.GPU_Veloc, octree.target_list,
                              octree.number_IP, octree.interaction_pairs, .05);

        }
#pragma omp single
        {
            phiTime = wcTime();
            UpSweep(octree.root);
            phiTime = wcTime() - phiTime;

            psiTime = wcTime();
            DownSweep(octree.root);
            psiTime = wcTime() - psiTime;
        }
#pragma omp barrier
#pragma omp single
        {
            gpuGetVelocities();
            dTime = wcTime() - dTime;

        }
#pragma omp barrier
        UpdateBodiesWithGPU(octree.numParticles, dt);

#pragma omp single
        {
            //ReSort(octree.root, octree.rootInfo);
            /*if(ReSort(octree.root, octree.rootInfo))
                RebuildTree(number_particles, precision, maximum_extent, minimum_extent);
                * */
        }
    }
    //if(dflag)
    //    printf("UP: %.3f, Down: %.3f, Direct: %.3f\n", phiTime, psiTime, dTime);
    fprintf(octree.output, "%.3f, %.3f, %.3f, %i", phiTime, psiTime, dTime, octree.numThreads);

}

/***********************************************************************
 *
 * set_affinity
 * @brief -  bind all threads in current parallel context to individual
 * cores openmp thread i will run on core i
 *
 **********************************************************************/
void set_affinity()
{
    int ret;
    cpu_set_t cpu_set;

    if (CPU_SETSIZE < omp_get_num_threads())
    {
        perror("Number of threads exceeed number of available cores");
        exit(1);
    }

    int myid = omp_get_thread_num(); /* omp thread id */
    int tid = syscall(__NR_gettid);  /* linux thread id */

#pragma omp barrier

    //printf("Setting CPU affinity on thread %d\n", myid);

    /* initialized bitmask */
    CPU_ZERO(&cpu_set);
    CPU_SET(myid, &cpu_set);

    ret = sched_setaffinity(tid, sizeof(cpu_set_t), &cpu_set);
    if (ret != 0)
    {
        perror("thread sched_setaffinity");
        exit(1);
    }

#pragma omp barrier

    return;
}


/***********************************************************************
 *
 * par_time
 * @brief - sum up whatever thread level vars we want
 *
 **********************************************************************/
void par_time(double *ccTot, double *phiTot, double *psiTot, double *diTot)
{
#pragma omp critical
    {
        //(void) fprintf(stderr,"Thread Level CreateCube() time: %8.3fs\n", cc_time);
        //(void) fprintf(stderr,"Thread Level ComputePhi() time: %8.3fs\n", phi_time);
        //(void) fprintf(stderr,"Thread Level ComputePsi() time: %8.3fs\n", psi_time);
        //(void) fprintf(stderr,"Thread Level direct time: %8.3fs\n", di_time);
        /* *ccTot += cc_time;
        *phiTot += phi_time;
        *psiTot += psi_time;
        *diTot += di_time;*/
        //(void) fprintf(stderr,"End Thread Info\n");
    }
}

/***********************************************************************
 *
 * AllPairs()
 * @brief - computes the forces on each particle via the direct, all
 * pairs method
 *
 **********************************************************************/
void AllPairs(unsigned long number_particles, double dt, double delta)
{
    int i, j;

#pragma omp parallel
    {
#pragma omp for private(i,j)
        for (i = 0; i < number_particles; i++)
        {
            for (j = 0; j < number_particles; j++)
            {
                ComputeVelocityDirect(i, j, delta);
            }
        }
//         UpdateBodies(number_particles, dt);
#pragma omp single
        {
            //ReSort(octree.root, octree.rootInfo);
        }
    }
}

/***********************************************************************
 *
 * UpdateBodies
 * @brief - updates position, velocity
 * soon will need dt component passed in
 *
 **********************************************************************/
void UpdateBodies(unsigned long number_particles, double dt)
{
    int i, idx;
#pragma omp for private(i)
    for (i = 0, idx = 0; i < number_particles; i++, idx += 3)
    {
        octree.bodies[i].position[0] += dt * octree.GPU_Veloc[idx];
        octree.bodies[i].position[1] += dt * octree.GPU_Veloc[idx+1];
        octree.bodies[i].position[2] += dt * octree.GPU_Veloc[idx+2];

        octree.bodies[i].position[0] = Modulus(octree.bodies[i].position[0], 256.0);
        octree.bodies[i].position[1] = Modulus(octree.bodies[i].position[1], 256.0);
        octree.bodies[i].position[2] = Modulus(octree.bodies[i].position[2], 256.0);

        //reset forces
        octree.GPU_Veloc[idx] = octree.GPU_Veloc[idx+1] = octree.GPU_Veloc[idx+2] = 0.0;

        /*assert((octree.bodies[i].x < 256.0)&&(octree.bodies[i].x >= 0) &&
                (octree.bodies[i].y < 256.0)&&(octree.bodies[i].y >= 0) &&
                (octree.bodies[i].z < 256.0)&&(octree.bodies[i].z >= 0));*/
    }
}

/***********************************************************************
 *
 * UpdateBodiesWithGPU
 * @brief - updates position, velocity
 *
 **********************************************************************/
void UpdateBodiesWithGPU(unsigned long number_particles, double dt)
{
    int i, idx;
#pragma omp for private(i)
    for (i = 0, idx = 0; i < number_particles; i++, idx += 3)
    {
        //combine CPU & GPU contributions
        octree.CPU_Veloc[idx] += octree.GPU_Veloc[idx];
        octree.CPU_Veloc[idx+1] += octree.GPU_Veloc[idx+1];
        octree.CPU_Veloc[idx+2] += octree.GPU_Veloc[idx+2];

        //update positions
        octree.bodies[i].position[0] += dt * octree.CPU_Veloc[idx];
        octree.bodies[i].position[1] += dt * octree.CPU_Veloc[idx+1];
        octree.bodies[i].position[2] += dt * octree.CPU_Veloc[idx+2];

        octree.bodies[i].position[0] = Modulus(octree.bodies[i].position[0], 256.0);
        octree.bodies[i].position[1] = Modulus(octree.bodies[i].position[1], 256.0);
        octree.bodies[i].position[2] = Modulus(octree.bodies[i].position[2], 256.0);

        //reset velocities
        octree.CPU_Veloc[idx] = octree.CPU_Veloc[idx+1] = octree.CPU_Veloc[idx+2] = 0.0;
        octree.GPU_Veloc[idx] = octree.GPU_Veloc[idx+1] = octree.GPU_Veloc[idx+2] = 0.0;

        if ((octree.bodies[i].position[0] >= 256.0) || (octree.bodies[i].position[0] < 0.0) ||
                (octree.bodies[i].position[1] >= 256.0) || (octree.bodies[i].position[1] < 0.0) ||
                (octree.bodies[i].position[2] >= 256.0) || (octree.bodies[i].position[2] < 0.0))
        {
            printf("x: %f, y:%f, z:%f\n", octree.bodies[i].position[0], octree.bodies[i].position[1], octree.bodies[i].position[2]);
        }
    }
}

/***********************************************************************
 *
 * PerformDirectInteractions(Node *node)
 * @brief Carry out all direct interactions in tree
 *
 **********************************************************************/
void PerformDirectInterations(Node *node,double delta)
{
    int id;

    if (!node->isParent)
    {
        CalculateDirectInteractions(node, node->list4, *node->list4Count,delta);
        PerformDirectWithinNode(node,delta);
        CalculateDirectInteractions(node, node->list1, *node->list1Count,delta);
        CalculateDirectInteractions(node, node->list3, *node->list3Count,delta);
    }

    if (node->isParent)
    {
        for (id = 0; id < MaxChildren; id++)
        {
            if (node->child[id]->pArrayLow != -1)
            {
#pragma omp task firstprivate(id)
                PerformDirectInterations(node->child[id],delta);
            }
        }
#pragma omp taskwait
    }
}

/***********************************************************************
 *
 * wcTime()
 * @brief - returns wall clock time
 *
 **********************************************************************/
double wcTime()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}

/***********************************************************************
 *
 * Plummer
 * @brief - generates plummer distribution as center of space in
 * particle array
 *
 **********************************************************************/
void Plummer(int N)
{
    int i, k;

    //dimenstions
    int ndim = 3;

    //used in distribution creation
    double mfrac = 0.999;
    double sqrt2 = 1.41421356237;
    double pi = 3.14159265359;
    double rsc = 3 * pi / 16;

    //Particle centerMass;
    //centerMass.position[0] = centerMass.position[1] = centerMass.position[2] = 0;

    double ri;

    for (i = 0; i < N; i++)
    {
        double temp = (((double) 1) / (pow((double)mfrac * randdouble(), ((double)1.5)))) - 1;
        ri = 1 / sqrt(temp);

        PickVector(rsc*ri, ndim, &octree.bodies[i].position[0], &octree.bodies[i].position[1], &octree.bodies[i].position[2]);

        octree.bodies[i].force[0] = .001;
        octree.bodies[i].force[1] = 0.0;
        octree.bodies[i].force[2] = 0.0;

        /*
        centerMass.position[0] += octree.bodies[i].position[0];
        centerMass.position[1] += octree.bodies[i].position[1];
        centerMass.position[2] += octree.bodies[i].position[2];
        */
    }
    /*
    centerMass.position[0] = centerMass.position[0]/N;
    centerMass.position[1] = centerMass.position[1]/N;
    centerMass.position[2] = centerMass.position[2]/N;
    */
    //shift to desired center
    for (i = 0; i < N; i++)
    {
        octree.bodies[i].position[0] += 128;
        octree.bodies[i].position[1] += 128;
        octree.bodies[i].position[2] += 128;
    }
}

/***********************************************************************
 *
 * PickVector is used to create a plummer distribution. See
 * void Plummer(int N, Particle * particles )
 *
 **********************************************************************/
void PickVector(double rad, int ndim, float *x, float *y, float *z)
{
    *x = 0;
    *y = 0;
    *z = 0;
    int tooBig = 1;
    double rsq;
    while (tooBig == 1)
    {
        rsq = 0;
        *x = 2 * randdouble() - 1;
        rsq = rsq + *x * *x;
        *y = 2 * randdouble() - 1;
        rsq = rsq + *y * *y;
        *z = 2 * randdouble() - 1;
        rsq = rsq + *z * *z;

        if (rsq > 1)
        {
            tooBig = 1;
        }
        else
        {
            tooBig = 0;
        }
    }

    *x = rad * *x / sqrt(rsq);
    *y = rad * *y / sqrt(rsq);
    *z = rad * *z / sqrt(rsq);
}

/***********************************************************************
 *
 * randDouble()
 * @brief - returns pseudo-random double
 *
 **********************************************************************/
double randdouble()
{
    return rand() / ((double)(RAND_MAX) + 1);
}

/***********************************************************************
 *
 * randFloat()
 * @brief - returns pseudo-random float
 *
 **********************************************************************/
float randfloat()
{
    return (float)rand() / ((float)(RAND_MAX) + 1);
}

/***********************************************************************
 *
 * PrintBodies
 *
 **********************************************************************/
void PrintBodies(int number_particles)
{
    int i;
    for (i = 0; i < number_particles; i++)
    {
        fprintf(stderr, "Particle %i, x: %.20f, y: %.20f, z: %.20f\n", i, octree.bodies[i].position[0], octree.bodies[i].position[1], octree.bodies[i].position[2]);
        //fprintf(stderr, "    fx %f, fy: %f, fz: %f\n", octree.CPU_Forces[i].x, octree.CPU_Forces[i].y, octree.CPU_Forces[i].z);
    }
}

/***********************************************************************
 *
 * Take a % b of floats, positive and neg
 *
 **********************************************************************/
float Modulus(float a, float b)
{
    if (a < 0)
        for (; a < 0; a += b);
    else
        for (; a >= b; a -= b);

    return a;
}


