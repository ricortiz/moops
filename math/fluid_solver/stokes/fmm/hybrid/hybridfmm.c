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
        for (id=0; id < 8; id++) {
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
                OuterToInner(node->list2[i], node, node->psi, 0);
                //number_interactions+=(node->pArrayHigh-node->pArrayLow+1)*(node->list2[i]->pArrayHigh-node->list2[i]->pArrayLow+1);
            }
        }
        DownShift(node, 0);
    }
    downTime+= wcTime() - stime;
    //recurse
    if (node->isParent)
    {
        for (id=0; id < 8; id++)
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
        for(j=0; j<8; j++){
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
        for(i = 0; i < 8; i++)
        {
            count = 0;
            if(octree.root->child[i]!=NULL)
            {
                colleaguesList = octree.root->child[i]->Colleagues;
                for(j=0; j<8; j++)
                {
                    if((octree.root->child[j]!=NULL)
                        &&(octree.root->child[j]!=octree.root->child[i])){
                        colleaguesList[count++] = octree.root->child[j];
                        }
                }
                for(k=0;k<8;k++){
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
                for(j=0;j<8;j++)
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
        for(i=0; i<8; i++)
        {
            if((node->parent->child[i]!=NULL)&&(node->parent->child[i]!=node)){
                colleaguesList[count++] = node->parent->child[i];
                assert(count<MaxColleagues);
            }
        }
        if((node->level<(octree.depth-1))&&(node->isParent)){
            for(k=0;k<8;k++){
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
        for(j=0; j<8; j++){
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
        for(i = 0; i<8; i++)
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
                for(j=0;j<8;j++){
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
        for(i = 0; i<8; i++)
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
        for(j = 0; j<8; j++){
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
                for(j = 0; j<8; j++){
                    CheckList3(node, node->Colleagues[i]->child[j]);
                }
            }
        }
    }
    
    if(node->isParent)//then has at least one child
    {
        for(i = 0; i<8; i++)
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
 * wcTime()
 * @brief - returns wall clock time
 *
 **********************************************************************/
double wcTime() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + 1E-6 * tv.tv_usec);
}
