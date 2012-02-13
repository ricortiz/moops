/*
 *  Include declarations.
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
 *  Global variables.
 */
FILE * gpuOut;
double upTime, downTime;
#pragma omp threadprivate(upTime, downTime)

void AssertList2(Node *node, Node* l2_node)
{
    float corners[24];
    
    double d;
    
    int i;
    
    float l2_r = octree.edge_length[l2_node->level]*0.5;
    
    corners[0] = l2_node->mid_x + l2_r;
    corners[1] = l2_node->mid_y + l2_r;
    corners[2] = l2_node->mid_z + l2_r;
    
    corners[3] = l2_node->mid_x + l2_r;
    corners[4] = l2_node->mid_y + l2_r;
    corners[5] = l2_node->mid_z - l2_r;
    
    corners[6] = l2_node->mid_x + l2_r;
    corners[7] = l2_node->mid_y - l2_r;
    corners[8] = l2_node->mid_z + l2_r;
    
    corners[9] = l2_node->mid_x - l2_r;
    corners[10] = l2_node->mid_y + l2_r;
    corners[11] = l2_node->mid_z + l2_r;
    
    corners[12] = l2_node->mid_x + l2_r;
    corners[13] = l2_node->mid_y - l2_r;
    corners[14] = l2_node->mid_z - l2_r;
    
    corners[15] = l2_node->mid_x - l2_r;
    corners[16] = l2_node->mid_y + l2_r;
    corners[17] = l2_node->mid_z - l2_r;
    
    corners[18] = l2_node->mid_x - l2_r;
    corners[19] = l2_node->mid_y - l2_r;
    corners[20] = l2_node->mid_z + l2_r;
    
    corners[21] = l2_node->mid_x - l2_r;
    corners[22] = l2_node->mid_y - l2_r;
    corners[23] = l2_node->mid_z - l2_r;
    
    float x,y,z;
    
    x = (corners[21]+corners[18])/2.0;
    y = (corners[22]+corners[19])/2.0;
    z = (corners[23]+corners[20])/2.0;
    
    d = sqrt((x-node->mid_x)*(x-node->mid_x) +
    (y-node->mid_y)*(y-node->mid_y)+
    (z-node->mid_z)*(z-node->mid_z)
    );
    if(d < (octree.edge_length[node->level]*0.5)*(sqrt(3)/2))
    {
        printf("Node: %.3f, %.3f, %.3f,%.3f \n             Bad list2: %.3f, %.3f, %.3f,%.3f\n",
               node->mid_x,node->mid_y,node->mid_z, octree.edge_length[node->level]*0.5,l2_node->mid_x,l2_node->mid_y,l2_node->mid_z, octree.edge_length[l2_node->level]*0.5 );
        exit(4);
    }
    for(i=0;i<8;i++)
    {
        d = sqrt((corners[3*i]-node->mid_x)*(corners[3*i]-node->mid_x) +
        (corners[3*i+1]-node->mid_y)*(corners[3*i+1]-node->mid_y)+
        (corners[3*i+2]-node->mid_z)*(corners[3*i+2]-node->mid_z)
        );
        //printf("D: %f\n",d);
        if(d < (octree.edge_length[node->level]*0.5)*(sqrt(3)/2))
        {
            printf("Node: %.3f, %.3f, %.3f,%.3f \n             Bad list2: %.3f, %.3f, %.3f,%.3f\n",
                   node->mid_x,node->mid_y,node->mid_z, octree.edge_length[node->level]*0.5,l2_node->mid_x,l2_node->mid_y,l2_node->mid_z, octree.edge_length[l2_node->level]*0.5 );
            exit(4);
        }
        
    }
}
/****************************************************************
 * 
 * returns binomial coefficient
 *
 ****************************************************************/
unsigned int Binomial(int n,int k)
{
    int i, j;
    unsigned int binomial;
    
    if (n < k)
    {
        return 0;
    }
    binomial=1;
    k=Min(k,n-k);
    j=n;
    for (i=1; i <= k; i++)
    {
        binomial=(unsigned long) ((binomial*(float) j/i)+0.5);
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
    
    if (node->isParent) {
        for (id=0; id < MaxChildren; id++) {
            if(node->child[id]->pArrayLow >= 0)
            {
                #pragma omp task firstprivate(id)
                UpSweep(node->child[id]);
            }
        }
        #pragma omp taskwait
        if(node->level > 0)
        {
            ShiftFromChildToParent(node);
        }
    }
    else
    {
        FormOuterExpansion(node);
    }
}

/***********************************************************************
 * 
 * DownSweep(Node *node)
 *
 **********************************************************************/
void DownSweep(Node *node)
{
    int id,i;
    double stime;
    
    stime = wcTime();
    if (!node->isParent)
    {
        //shift in list 2
        for(i = 0; i < *(node->list2Count); i++)
        {
            if(node->list2[i]->pArrayLow>=0){
                //Direct(node, node->list2[i]);
                OuterToInner(node->list2[i], node, node->psi, True);
                //number_interactions+=(node->pArrayHigh-node->pArrayLow+1)*(node->list2[i]->pArrayHigh-node->list2[i]->pArrayLow+1);
            }
        }
        DownShift(node, True);
        for(i=0;i<*node->list3Count;i++)
        {
            if((node->list3[i]->isParent)&&(node->list3[i]->pArrayLow>=0)){
                //Direct(node, node->list3[i]);
                OuterToInner(node->list3[i], node, node->psi, True);
                //number_interactions+=(node->pArrayHigh-node->pArrayLow+1)*(node->list3[i]->pArrayHigh-node->list3[i]->pArrayLow+1);
            }
        }
        ApplyLocalExpansion(node);
    }
    else if (node->level > 1)
    {
        for(i = 0; i < *(node->list2Count); i++)
        {
            if(node->list2[i]->pArrayLow>=0){
                //Direct(node, node->list2[i]);
                OuterToInner(node->list2[i], node, node->psi, False);
                //number_interactions+=(node->pArrayHigh-node->pArrayLow+1)*(node->list2[i]->pArrayHigh-node->list2[i]->pArrayLow+1);
            }
        }
        DownShift(node, False);
    }
    downTime+= wcTime() - stime;
    //recurse
    if (node->isParent)
    {
        for (id=0; id < MaxChildren; id++)
        {
            if(node->child[id]->pArrayLow!=-1)
            {
                #pragma omp task firstprivate(id)
                DownSweep(node->child[id]);
            }
        }
        #pragma omp taskwait
    }
    
}

void CheckColleagues(Node *node, Node *candidate){
    int j;
    if((candidate->level < node->level) && candidate->isParent){
        for(j=0; j<MaxChildren; j++){
            if((candidate->child[j]!=NULL)&&(IsBorder(node, candidate->child[j])==1)){
                CheckColleagues(node, candidate->child[j]);
            }
        }
    }
    else{
        if(IsBorder(node, candidate)==1){
            node->Colleagues[*node->colleagueCount] = candidate;
            *node->colleagueCount = *node->colleagueCount + 1;
        }
    }
}
/**
 * void BuildColleagues(Node *node){
 *    int i,j,k,count;
 *    Node **colleaguesList, **srcColleaguesList;
 *    if(node == octree.root)
 *    {
 *        for(i = 0; i < MaxChildren; i++)
 *        {
 *            count = 0;
 *            if(octree.root->child[i]!=NULL)
 *            {
 *                colleaguesList = octree.root->child[i]->Colleagues;
 *                for(j=0; j<MaxChildren; j++)
 *                {
 *                    if((octree.root->child[j]!=NULL)
 *                    &&(octree.root->child[j]!=octree.root->child[i])){
 *                        colleaguesList[count++] = octree.root->child[j];
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
void BuildColleagues(Node * node){
    int i, j, k, count;
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
        Node *stack[MaxList1Stack], *temp;
        int top;
        top = -1;
        for(i = 0; i<(MaxList1Stack); i++)
        {
            stack[i] = NULL;
        }
        
        srcColleaguesList = node->parent->Colleagues;
        colleaguesList = node->Colleagues;
        count = 0;
        for(i = 0; i < (MaxList3); i++) //check parents colleagues
        {
            if(srcColleaguesList[i]!=NULL)
            {
                //put any node that borders @node on the stack
                //before adding to node->BorderList we will check to
                //see if it can be subdivided more
                if(IsBorder(node, srcColleaguesList[i]) == 1)
                {
                    stack[++top] = srcColleaguesList[i];
                    assert(top<MaxList1Stack);
                }
            }
        }
        while(top!=-1)
        {
            //if it has children and is at earlier level than node
            if(stack[top]->isParent&&(stack[top]->level<node->level))
            {
                //take off stack since only interested in its children now
                temp = stack[top--];
                for(j=0;j<MaxChildren;j++)
                {
                    if(temp->child[j]!=NULL)
                    {
                        if(IsBorder(node,temp->child[j])==1)//only add if borders
                        {
                            stack[++top] = temp->child[j];
                            assert(top<MaxList1Stack);
                        }
                    }
                }
            }
            else{ //put in colleagues list
                //assert(stack[top]->level<=node->level);
                colleaguesList[count++] = stack[top];
                assert(count<MaxColleagues);
                stack[top--] = NULL;
            }
            
        }
        for(i=0; i<MaxChildren; i++)
        {
            if((node->parent->child[i]!=NULL)&&(node->parent->child[i]!=node)){
                colleaguesList[count++] = node->parent->child[i];
                assert(count<MaxColleagues);
            }
        }
        if((node->level<(octree.depth-1))&&(node->isParent)){
            for(k=0;k<MaxChildren;k++){
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
void CheckList1(Node *node, Node *candidate){
    int j;
    if(candidate->isParent){
        for(j=0; j<MaxChildren; j++){
            if(IsBorder(node, candidate->child[j])==1){
                CheckList1(node, candidate->child[j]);
            }
        }
    }
    else{
        node->list1[*node->list1Count] = candidate;
        *node->list1Count = *node->list1Count + 1;
        assert((*node->list1Count)<MaxList1);
    }
}

/*********************************************************************
 * 
 * Builds List1 for a node and its children
 *
 ********************************************************************/
void BuildList1(Node *node){
    int i;
    if(!node->isParent)
    {
        //check each colleague
        for(i = 0; i<MaxColleagues; i++)
        {
            //then seen all colleagues
            if(node->Colleagues[i]==NULL)
            {
                break;
            }
            else{
                CheckList1(node, node->Colleagues[i]);
            }
        }
    }
    
    if(node->isParent)
    {
        for(i = 0; i<MaxChildren; i++)
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
void BuildList2(Node *node){
    int i,j;
    
    for(i=0;i<MaxColleagues;i++){
        if(node->parent->Colleagues[i]==NULL){
            break;
        }
        else if(node->parent->level == node->parent->Colleagues[i]->level){ //then true colleague
            if(node->parent->Colleagues[i]->isParent){
                for(j=0;j<MaxChildren;j++){
                    if(IsBorder(node,node->parent->Colleagues[i]->child[j])!=1)
                    {
                        node->list2[*node->list2Count] = node->parent->Colleagues[i]->child[j];
                        *node->list2Count = *node->list2Count + 1;
                        assert((*node->list2Count)<MaxList3);
                    }
                }
            }
        }
    }
    //descend
    if(node->isParent)
    {
        for(i = 0; i<MaxChildren; i++)
        {
            #pragma omp task firstprivate(i)
            BuildList2(node->child[i]);
        }
        #pragma omp taskwait
    }
}

void CheckList3(Node * node, Node * candidate){
    int j;
    if(IsBorder(node, candidate)!=1) //if not adjacent
    {
        //put in list3 and by result list4
        node->list3[*node->list3Count] = candidate;
        *node->list3Count = *node->list3Count + 1;
        assert((*node->list3Count)<MaxList3);
        candidate->list4[*candidate->list4Count] = node;
        *candidate->list4Count = *candidate->list4Count + 1;
        assert((*node->list4Count)<MaxList4);
    }
    else if(candidate->isParent){
        for(j = 0; j<MaxChildren; j++){
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
    int i,j;
    
    if(!node->isParent)
    {
        //for each colleague
        for(i = 0; i<MaxColleagues; i++)
        {
            //then seen all colleagues
            if(node->Colleagues[i]==NULL)
            {
                break;
            }
            else if(node->Colleagues[i]->isParent){ //each child of colleague is potential list3 node
                for(j = 0; j<MaxChildren; j++){
                    CheckList3(node, node->Colleagues[i]->child[j]);
                }
            }
        }
    }
    
    if(node->isParent)//then has at least one child
    {
        for(i = 0; i<MaxChildren; i++)
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
 **********************************************************************
 int IsBorder(Node * constructListFor, Node * toCheck){
     
     float toCheckMidX, toCheckMidY, toCheckMidZ, constructMidX, constructMidY, constructMidZ,
     minusCheckX, plusCheckX, minusCheckY, plusCheckY, minusCheckZ, plusCheckZ, minusConstX,
     plusConstX, minusConstY, plusConstY, minusConstZ, plusConstZ;
     double toCheckBisect, constructBisect, eps;
     
     int toReturn;
     
     eps = 1e-10;
     
     toCheckMidX = toCheck->mid_x;
     toCheckMidY = toCheck->mid_y;
     toCheckMidZ = toCheck->mid_z;
     toCheckBisect = octree.edge_length[toCheck->level]*0.5;
     
     minusCheckX = toCheckMidX - toCheckBisect;
     plusCheckX = toCheckMidX + toCheckBisect;
     minusCheckY = toCheckMidY - toCheckBisect;
     plusCheckY = toCheckMidY + toCheckBisect;
     minusCheckZ = toCheckMidZ - toCheckBisect;
     plusCheckZ = toCheckMidZ + toCheckBisect;
     
     constructMidX = constructListFor->mid_x;
     constructMidY = constructListFor->mid_y;
     constructMidZ = constructListFor->mid_z;
     constructBisect = octree.edge_length[constructListFor->level]*0.5;
     
     minusConstX = constructMidX - constructBisect;
     plusConstX = constructMidX + constructBisect;
     minusConstY = constructMidY - constructBisect;
     plusConstY = constructMidY + constructBisect;
     minusConstZ = constructMidZ - constructBisect;
     plusConstZ = constructMidZ + constructBisect;
     
     if((AbsoluteValue(minusConstX - plusCheckX)<eps)||(AbsoluteValue(plusConstX - minusCheckX)<eps))
     {
         if( ((plusConstY>minusCheckY)||(AbsoluteValue(plusConstY-minusCheckY)<eps))
             &&((minusConstY<plusCheckY)||(AbsoluteValue(minusConstY-plusCheckY)<eps))
             &&((minusConstZ<plusCheckZ)||(AbsoluteValue(minusConstZ-plusCheckZ)<eps))
             &&((plusConstZ>minusCheckZ)||(AbsoluteValue(plusConstZ-minusCheckZ)<eps))
             )
             {
                 return 1;
             }
     }
     else if((AbsoluteValue(minusConstZ-plusCheckZ)<eps)||(AbsoluteValue(plusConstZ-minusCheckZ)<eps))
     {
         if
             ( ((plusConstX>minusCheckX)||(AbsoluteValue(plusConstX-minusCheckX)<eps))
             &&((minusConstX<plusCheckX)||(AbsoluteValue(minusConstX-plusCheckX)<eps))
             &&((minusConstY<plusCheckY)||(AbsoluteValue(minusConstY-plusCheckY)<eps))
             &&((plusConstY>minusCheckY)||(AbsoluteValue(plusConstY-minusCheckY)<eps))
             )
             {
                 return 1;
             }
     }
     else if((AbsoluteValue(minusConstY-plusCheckY)<eps)||(AbsoluteValue(plusConstY-minusCheckY)<eps))
     {
         if( ((plusConstX>minusCheckX)||(AbsoluteValue(plusConstX-minusCheckX)<eps))
             &&((minusConstX<plusCheckX)||(AbsoluteValue(minusConstX-plusCheckX)<eps))
             &&((plusConstZ>minusCheckZ)||(AbsoluteValue(plusConstZ-minusCheckZ)<eps))
             &&((minusConstZ<plusCheckZ)||(AbsoluteValue(minusConstZ-plusCheckZ)<eps))
             )
             {
                 return 1;
             }
     }
     return 0;
 }*/

/***********************************************************************
 * 
 * IsBorder(Node * constructListFor, Node * toCheck) checks if toCheck
 * borders constructListFor and returns 1 if so, 0 otherwise.
 *
 **********************************************************************
 int IsBorder(Node * constructListFor, Node * toCheck){
     
     float toCheckMidX, toCheckMidY, toCheckMidZ, constructMidX, constructMidY, constructMidZ,
     minusCheckX, plusCheckX, minusCheckY, plusCheckY, minusCheckZ, plusCheckZ, minusConstX,
     plusConstX, minusConstY, plusConstY, minusConstZ, plusConstZ;
     double toCheckBisect, constructBisect;
     
     int toReturn;
     
     toCheckMidX = toCheck->mid_x;
     toCheckMidY = toCheck->mid_y;
     toCheckMidZ = toCheck->mid_z;
     toCheckBisect = octree.edge_length[toCheck->level]*0.5;
     
     minusCheckX = toCheckMidX - toCheckBisect;
     plusCheckX = toCheckMidX + toCheckBisect;
     minusCheckY = toCheckMidY - toCheckBisect;
     plusCheckY = toCheckMidY + toCheckBisect;
     minusCheckZ = toCheckMidZ - toCheckBisect;
     plusCheckZ = toCheckMidZ + toCheckBisect;
     
     constructMidX = constructListFor->mid_x;
     constructMidY = constructListFor->mid_y;
     constructMidZ = constructListFor->mid_z;
     constructBisect = octree.edge_length[constructListFor->level]*0.5;
     
     minusConstX = constructMidX - constructBisect;
     plusConstX = constructMidX + constructBisect;
     minusConstY = constructMidY - constructBisect;
     plusConstY = constructMidY + constructBisect;
     minusConstZ = constructMidZ - constructBisect;
     plusConstZ = constructMidZ + constructBisect;
     
     if((minusConstX == plusCheckX)||(plusConstX == minusCheckX))
     {
         if((plusConstY>=minusCheckY)&&(minusConstY<=plusCheckY)
             &&(minusConstZ<=plusCheckZ)&&(plusConstZ>=minusCheckZ))
             {
                 return 1;
             }
     }
     else if((minusConstZ==plusCheckZ)||(plusConstZ==minusCheckZ))
     {
         if((plusConstX>=minusCheckX)&&(minusConstX<=plusCheckX)
             &&(minusConstY<=plusCheckY)&&(plusConstY>=minusCheckY))
             {
                 return 1;
             }
     }
     else if((minusConstY==plusCheckY)||(plusConstY==minusCheckY))
     {
         if((plusConstX>=minusCheckX)&&(minusConstX<=plusCheckX)
             &&(plusConstZ>=minusCheckZ)&&(minusConstZ<=plusCheckZ))
             {
                 return 1;
             }
     }
     return 0;
 }
 */
int IsBorder(Node *node1, Node *node2)
{
    float extents[2] = {octree.edge_length[node1->level]*.5,octree.edge_length[node2->level]*.5};
    float point[4][3] =
    {
        {node1->mid_x - extents[0],node1->mid_y - extents[0],node1->mid_z - extents[0]},
        {node1->mid_x + extents[0],node1->mid_y + extents[0],node1->mid_z + extents[0]},
        {node2->mid_x - extents[1],node2->mid_y - extents[1],node2->mid_z - extents[1]},
        {node2->mid_x + extents[1],node2->mid_y + extents[1],node2->mid_z + extents[1]}
    };
    
    int result = ((point[0][0] - point[2][0] <= 1e-6 && point[2][0] - point[1][0] <= 1e-6) || ( point[0][0] - point[3][0] <= 1e-6 && point[3][0] - point[1][0] <= 1e-6));
    result = result && (( point[0][1] - point[2][1] <= 1e-6 && point[2][1] - point[1][1] <= 1e-6) || ( point[0][1] - point[3][1] <= 1e-6 && point[3][1] - point[1][1] <= 1e-6));
    result = result && (( point[0][2] - point[2][2] <= 1e-6 && point[2][2] - point[1][2] <= 1e-6) || ( point[0][2] - point[3][2] <= 1e-6 && point[3][2] - point[1][2] <= 1e-6));
    
    return result;
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
    
    if (CPU_SETSIZE < omp_get_num_threads()) {
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
    if (ret != 0) {
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
void par_time(double *ccTot, double *phiTot, double *psiTot, double *diTot){
    #pragma omp critical
    {
        //(void) fprintf(stderr,"Thread Level CreateCube() time: %8.3fs\n", cc_time);
        //(void) fprintf(stderr,"Thread Level ComputePhi() time: %8.3fs\n", phi_time);
        //(void) fprintf(stderr,"Thread Level ComputePsi() time: %8.3fs\n", psi_time);
        //(void) fprintf(stderr,"Thread Level direct time: %8.3fs\n", di_time);
        /* *ccTot += cc_time;
         * phiTot += phi_time;
         *psiTot += psi_time;
         *diTot += di_time;*/
        //(void) fprintf(stderr,"End Thread Info\n");
    }
}




/***********************************************************************
 * 
 * wcTime()
 * @brief - returns wall clock time
 *
 **********************************************************************/
double wcTime() {
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
void Plummer(int N){
    int i,k;
    
    //dimenstions
    int ndim = 3;
    
    //used in distribution creation
    double mfrac = 0.999;
    double sqrt2 = 1.41421356237;
    double pi = 3.14159265359;
    double rsc = 3 * pi/16;
    
    //Particle centerMass;
    //centerMass.position[0] = centerMass.position[1] = centerMass.position[2] = 0;
    
    double ri;
    
    for(i=0; i<N; i++)
    {
        double temp = (((double) 1 )/(pow((double)mfrac*randdouble(), ((double)1.5)))) - 1;
        ri = 1/sqrt(temp);
        
        PickVector(rsc*ri, ndim, &octree.bodies[i].position[0],&octree.bodies[i].position[1],&octree.bodies[i].position[2]);
        
        octree.bodies[i].force[0] = .001;
        octree.bodies[i].force[1] = 0.0;
        octree.bodies[i].force[2] = 0.0;
        
        /*
         *        centerMass.position[0] += octree.bodies[i].position[0];
         *        centerMass.position[1] += octree.bodies[i].position[1];
         *        centerMass.position[2] += octree.bodies[i].position[2];
         */
    }
    /*
     *    centerMass.position[0] = centerMass.position[0]/N;
     *    centerMass.position[1] = centerMass.position[1]/N;
     *    centerMass.position[2] = centerMass.position[2]/N;
     */
    //shift to desired center
    for(i=0; i<N; i++)
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
        *x = 2*randdouble() - 1;
        rsq = rsq + *x * *x;
        *y = 2*randdouble() - 1;
        rsq = rsq + *y * *y;
        *z = 2*randdouble() - 1;
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
    
    *x = rad * *x/sqrt(rsq);
    *y = rad * *y/sqrt(rsq);
    *z = rad * *z/sqrt(rsq);
}

/***********************************************************************
 * 
 * randDouble()
 * @brief - returns pseudo-random double
 *
 **********************************************************************/
double randdouble()
{
    return rand()/((double)(RAND_MAX)+1);
}

/***********************************************************************
 * 
 * randFloat()
 * @brief - returns pseudo-random float
 *
 **********************************************************************/
float randfloat()
{
    return (float)rand()/((float)(RAND_MAX)+1);
}

/***********************************************************************
 * 
 * PrintBodies
 *
 **********************************************************************/
void PrintBodies(int number_particles)
{
    int i;
    for (i=0; i<number_particles; i++)
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
    if(a<0)
        for(; a < 0; a += b);
        else
            for(; a >= b; a -= b);
            
            return a;
}
