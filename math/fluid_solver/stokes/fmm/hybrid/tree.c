#include "hybridfmm.h"
#include <assert.h>

int nextInfo,threadId, maxDepth;
NodeInfo * info_array, *root_info;
Node **DI_nodes;
unsigned long nextDI_node=0;
unsigned int *tempTargetList, *tempNum_IP;

#pragma omp threadprivate (info_array, nextInfo, maxDepth, threadId, DI_nodes, nextDI_node)

/***********************************************************************
 *
 * CreateOctree()
 * @brief - build the octree and the corresponding lists
 *
 **********************************************************************/
void CreateOctree(unsigned long number_particles, int precision, double maximum_extent, double minimum_extent)
{    
    int i, failure = True;
    double startTime;
    
    treeVarSetup(number_particles, precision, maximum_extent, minimum_extent);    
    octree.nodeInfoCoeff = 5;
    octree.treeRebuilds = -1;    
    
    //one-time setup
    #pragma omp parallel
    {
        #pragma omp single
        {
            octree.nodeInfoDistribution = number_particles/omp_get_num_threads();
            octree.threadFailures = (int *) malloc(omp_get_num_threads()*sizeof(int));
            octree.threads = omp_get_num_threads();
        }
    }
    
    while(failure){
        octree.treeRebuilds++;
        octree.nodeInfoCoeff = octree.nodeInfoCoeff + 1;
        constructTreeInfo(number_particles,minimum_extent);
        failure = False;
        for(i=0; i<octree.threads; i++)
        {
            failure = failure || octree.threadFailures[i];                
        }
        if(failure)
        {
            #pragma omp parallel
            {
                free((NodeInfo *) info_array);
                info_array = (NodeInfo *) NULL;
            }
        }
    }
 
    octree.depth = 0;
    #pragma omp parallel
    {
        #pragma omp critical
        {
            octree.depth = (maxDepth > octree.depth) ? maxDepth:octree.depth;
        }
        #pragma omp barrier
        #pragma omp single
        {
            fprintf(stderr, "Total Nodes: %i\n", octree.rootInfo->subtree_size);
            octree.nodes = (Node *) malloc(octree.rootInfo->subtree_size*sizeof(Node));
            for(i=0;i<4;i++)
            {
                octree.PHIs[i] = (float *) malloc(octree.rootInfo->subtree_size*octree.coefficients*sizeof(float));
                octree.PSIs[i] = (float *) malloc(octree.rootInfo->subtree_size*octree.coefficients*sizeof(float));
            } 

            octree.counts = (int *) malloc((octree.rootInfo->subtree_size*5)*sizeof(int));
            if((octree.nodes == NULL)||(octree.PHIs[0] == NULL) ||(octree.PHIs[1] == NULL)||(octree.PHIs[2] == NULL)||(octree.PHIs[3] == NULL)
            || (octree.PSIs[0] == NULL) || (octree.PSIs[1] == NULL)|| (octree.PSIs[2] == NULL)|| (octree.PSIs[3] == NULL) ||(octree.counts == NULL))
            {
                printf("Not enough mem for nodes/phis/psis\n");
                exit(4);
            }

            for(i=0;i<4;i++)
            {
                memset(octree.PHIs[i], 0.0, octree.rootInfo->subtree_size*octree.coefficients*sizeof(float));
                memset(octree.PSIs[i], 0.0, octree.rootInfo->subtree_size*octree.coefficients*sizeof(float));
            } 

            memset(octree.counts, 0, octree.rootInfo->subtree_size*5*sizeof(int));

            /*
             * have to initialize root prior to CreateSubtree to satisfy the
             * invariant that the parent it is called on has already been init
            */
            octree.root = InitializeNode(0, 0, NULL, 
                octree.rootInfo->mid_x, octree.rootInfo->mid_y, octree.rootInfo->mid_z, 
                octree.rootInfo->low, octree.rootInfo->high, 0, octree.rootInfo->isParent);
            octree.root->parent = octree.root;
            CreateSubtree(octree.root, octree.rootInfo, 1);

            startTime = wcTime();
            //build lists
            BuildColleagues(octree.root);  
            BuildList1(octree.root);
            BuildList2(octree.root);
            BuildList3And4(octree.root);
            octree.numThreads = omp_get_num_threads();
        }
    }

    /*
     * same structure for creating octree.DInodes array as for
     * building tree
     */
    octree.DI_coeff = 6.0;
    octree.treeRebuilds = -1;
    failure = True;

    startTime = wcTime();
    while(failure){
        octree.treeRebuilds++;
        octree.DI_coeff*= 2.0;
        constructDINodes(number_particles,(double)octree.rootInfo->subtree_size);
        
        failure = False;
        for(i=0; i<octree.threads; i++)
        {
            failure = failure || octree.threadFailures[i];                
        }
        if(failure)
        {
            #pragma omp parallel
            {
                free((Node **) DI_nodes);
                DI_nodes = (Node **) NULL;
            }
        }
    }
    //printf("Rebuilds: %i\n",  octree.treeRebuilds);
    BuildGPUArrays();

    printf("Tree depth: %lu\n", octree.depth);
}

void FreeOctree()
{
    int i;
    
    //tree and nodes
    free((Node *) octree.nodes);
    for(i=0;i<4;i++)
    {
        free((float *) octree.PHIs[i]);
        free((float *) octree.PSIs[i]);
        octree.PHIs[i] = octree.PSIs[i] = (float *) NULL;
    }
    free((int *) octree.counts);
    free((float *) octree.ijk_factorial);
    free((float *) octree.ijk_binomial);
    free((int *) octree.threadFailures);
    
    octree.nodes = (Node *) NULL;
    octree.ijk_factorial = octree.ijk_binomial = (float *) NULL;
    octree.counts = octree.threadFailures = (int *) NULL;
    
    //GPU data
    free((Node **) octree.DInodes);
    free((Node **) octree.leaf_DI_nodes);
    free((unsigned int *) tempTargetList);
    free((unsigned int *) tempNum_IP);
    free((unsigned int *) octree.interaction_pairs);
    
    octree.DInodes = octree.leaf_DI_nodes = (Node **) NULL;
    octree.target_list = octree.number_IP = octree.interaction_pairs = (unsigned int *) NULL;
}
/***********************************************************************
 *
 * constructDINodes()
 * @brief - sets up thread level DI_node arrays, then calls
 *  ConstructDINodesRec() (Rec for recursive) to fill up these thread
 *  level arrays. Then dumps all these arrays into the one main array at
 *  octree.DInodes
 *
 **********************************************************************/
void constructDINodes(unsigned long number_particles, double numNodes)
{    
    octree.DI_nodeTCount = ceil((octree.DI_coeff*numNodes)/(double)octree.numThreads);
    #pragma omp parallel
    {
        DI_nodes = (Node **) malloc(octree.DI_nodeTCount*sizeof(Node *));
        if(DI_nodes == NULL)
        {
            printf("Ran out of memory on threadlevel di node array.\n");
            exit(4);
        }
        nextDI_node = 0;
        octree.threadFailures[threadId] = False;
        #pragma omp barrier
        #pragma omp single
        {
            ConstructDINodesRec(octree.root);
        }
    }
}

/***********************************************************************
 * 
 * PushDownDirectWork
 * 
 * To avoid race conditions on GPU, if a parent does direct work, push it
 * down to any children, grandchildren, etc so that only leaf nodes do 
 * direct work
 * 
 **********************************************************************/
 void PushDownDirectWork(Node * parent)
 {
    int i, id;
    
    //assert(parent->isParent);
    
    for(i=0;i<*parent->list4Count;i++)
    {
        //if > maxBodies then done as expansion on CPU
        if((!parent->list4[i]->isParent)&&(parent->list4[i]->pArrayLow!=-1))
        {
            for(id=0;id<MaxChildren;id++)
            {
                if(parent->child[id]->pArrayLow!=-1)
                {
                    //printf("before %i\n", *parent->child[id]->list4Count);
                    parent->child[id]->list4[*parent->child[id]->list4Count] = parent->list4[i];
                    *parent->child[id]->list4Count = *parent->child[id]->list4Count + 1;
                    //printf("after %i\n", *parent->child[id]->list4Count);
                }
            }
        }
    }
    
    for(id=0;id<MaxChildren;id++)
    {
        //if child is parent and not a DI node, then continue to push
        if((parent->child[id]->isParent)&&(parent->child[id]->DI_tag==-1))
        {
            PushDownDirectWork(parent->child[id]);
        }
    }
 }

/***********************************************************************
 *
 *BuildGPUArrays()
 * @brief - constructs arrays to be handed off to GPU for direct
 * interaction computation
 *
 **********************************************************************/
void BuildGPUArrays()
{
    unsigned int DIN_prefix[octree.numThreads+1], i,j, *preSumNum_IP;
    double startTime;
    
    DIN_prefix[0] = 0;
    
    unsigned int numleafDISeen=0;
    
    #pragma omp parallel firstprivate(i) reduction(+ : numleafDISeen)
    {
        DIN_prefix[threadId+1] = nextDI_node;
        #pragma omp barrier
        #pragma omp single
        { 
            for(i=0;i<octree.numThreads;i++)
            {
                DIN_prefix[i+1] += DIN_prefix[i];
            }
            
            octree.number_DInodes = DIN_prefix[octree.numThreads];

            octree.DInodes = (Node **) malloc(octree.number_DInodes*sizeof(Node *));
            if(octree.DInodes== NULL)
            {
                printf("Not enough mem for DI nodes array\n");
                exit(4);
            }        
        }
        #pragma omp barrier
        for(i=0;i<nextDI_node;i++)
        {
            octree.DInodes[i+DIN_prefix[threadId]] = DI_nodes[i];
            if(!DI_nodes[i]->isParent)
                numleafDISeen++;
            //marker to show it is DI node for now
            octree.DInodes[i+DIN_prefix[threadId]]->DI_tag = -2;
            //printf("Low: %i, High: %i, isParent: %i\n", octree.DInodes[i+DIN_prefix[threadId]]->pArrayLow, octree.DInodes[i+DIN_prefix[threadId]]->pArrayHigh, octree.DInodes[i+DIN_prefix[threadId]]->isParent);
        }
    }
    
    octree.numLeafDInodes = numleafDISeen;
    octree.leaf_DI_nodes = (Node **) malloc(octree.numLeafDInodes*sizeof(Node *));
    if(octree.leaf_DI_nodes== NULL)
    {
        printf("Not enough mem for leaf DI nodes array\n");
        exit(4);
    }

    int idx = 0;
    for(j=0; j<octree.number_DInodes; j++)
    {
        if(octree.DInodes[j]->isParent)
        {
            PushDownDirectWork(octree.DInodes[j]);
        }
        else
        {
            octree.leaf_DI_nodes[idx] = octree.DInodes[j];
            octree.leaf_DI_nodes[idx]->DI_tag = idx++;
        }
    }
    
    
    tempTargetList = (unsigned int *) malloc(((2*octree.numLeafDInodes)*sizeof(unsigned int*)) + 2*MEMORY_ALIGNMENT);
    if(tempTargetList== NULL)
    {
        printf("Not enough mem for target list\n");
        exit(4);
    }
    octree.target_list = (unsigned int *) ALIGN_UP( tempTargetList, MEMORY_ALIGNMENT );
    
    tempNum_IP = (unsigned int *) malloc((octree.numLeafDInodes*sizeof(unsigned int)) + 2*MEMORY_ALIGNMENT);
    if(tempNum_IP == NULL)
    {
        printf("Not enough mem for number IP array\n");
        exit(4);
    }
    octree.number_IP = (unsigned int *) ALIGN_UP( tempNum_IP, MEMORY_ALIGNMENT );

    preSumNum_IP = (unsigned int *) malloc(octree.numLeafDInodes*sizeof(unsigned int));
    if(preSumNum_IP == NULL)
    {
        printf("Not enough mem for prefix sum of number IP array\n");
        exit(4);
    }

    for(j=0; j<octree.numLeafDInodes; j++)
    {
        octree.target_list[2*j] = octree.leaf_DI_nodes[j]->pArrayLow;
        octree.target_list[2*j+1] = octree.leaf_DI_nodes[j]->pArrayHigh - octree.leaf_DI_nodes[j]->pArrayLow +1;
    }

    for(j=0; j<octree.numLeafDInodes; j++)
    {
        octree.number_IP[j] = 1; //interact with self
        
        for(i=0; i<*octree.leaf_DI_nodes[j]->list1Count; i++)
        {
            if(octree.leaf_DI_nodes[j]->list1[i]->pArrayLow != -1)
            {
                octree.number_IP[j]+=1;
            }
        }
        for(i=0; i<*octree.leaf_DI_nodes[j]->list3Count; i++)
        {
            if((!octree.leaf_DI_nodes[j]->list3[i]->isParent) && (octree.leaf_DI_nodes[j]->list3[i]->pArrayLow != -1))
            {
                octree.number_IP[j]+=1;
            }
        }
        for(i=0; i<*octree.leaf_DI_nodes[j]->list4Count; i++)
        {
            if(octree.leaf_DI_nodes[j]->list4[i]->pArrayLow != -1)
            {
                octree.number_IP[j]+=1;
            }
        }
    }

    preSumNum_IP[0] = octree.number_IP[0];
    for(j=1; j<octree.numLeafDInodes; j++)
    {
        preSumNum_IP[j] = preSumNum_IP[j-1] + octree.number_IP[j];
    }
    octree.total_interaction_pairs = preSumNum_IP[octree.numLeafDInodes-1];  
    octree.interaction_pairs = (unsigned int *) malloc(octree.total_interaction_pairs*sizeof(unsigned int));
    if(octree.interaction_pairs== NULL)
    {
        printf("Not enough mem for the huge interaction pairs array\n");
        exit(4);
    }

    int currIdx = 0;
    for(j=0;j<octree.numLeafDInodes;j++)
    {        
        octree.interaction_pairs[currIdx++] = octree.leaf_DI_nodes[j]->DI_tag;

        for(i=0; i<*octree.leaf_DI_nodes[j]->list1Count; i++)
        {
            if(octree.leaf_DI_nodes[j]->list1[i]->pArrayLow != -1)
            {
                octree.interaction_pairs[currIdx++] = octree.leaf_DI_nodes[j]->list1[i]->DI_tag;
            }
        }

        for(i=0; i<*octree.leaf_DI_nodes[j]->list3Count; i++)
        {
            if((!octree.leaf_DI_nodes[j]->list3[i]->isParent) && (octree.leaf_DI_nodes[j]->list3[i]->pArrayLow != -1))
            {
                octree.interaction_pairs[currIdx++] = octree.leaf_DI_nodes[j]->list3[i]->DI_tag;
            }
        }

        for(i=0; i<*octree.leaf_DI_nodes[j]->list4Count; i++)
        {
            if(octree.leaf_DI_nodes[j]->list4[i]->pArrayLow != -1)
            {
                octree.interaction_pairs[currIdx++] = octree.leaf_DI_nodes[j]->list4[i]->DI_tag;
            }
        }
    }
    
    free((unsigned int *) preSumNum_IP);
    preSumNum_IP = (unsigned int *) NULL;
}
/***********************************************************************
 *
 * ConstructDINodesRec()
 * @brief - moves through the tree adding parents and leaves to thread
 *  level DI_nodes array if applicable and then task recursing
 *
 **********************************************************************/
void ConstructDINodesRec(Node *node)
{
    int id;

    if(nextDI_node>=octree.DI_nodeTCount)
    {
        octree.threadFailures[threadId] = True;
    }
    else if(!octree.threadFailures[threadId])
    {
        if(node->isParent)
        {
            if((*node->list4Count)>0)
            {
                DI_nodes[nextDI_node++] = node;
            }
            for(id=0;id<MaxChildren;id++)
            {
                if(node->child[id]->pArrayLow!=-1)
                {
                    #pragma omp task firstprivate(id)
                    ConstructDINodesRec(node->child[id]);
                }
            }
            #pragma omp taskwait
        }
        else{
            //assert(node->pArrayLow!=-1);
            DI_nodes[nextDI_node++] = node;
        }
    }
}
/***********************************************************************
 * 
 * Given an already initialized parent, build the subtree for each 
 * child of that parent (if there are any) 
 * 
 * Builds depth first
 * 
 * parent: the root of subtree to be created
 * 
 * parentInfo: NodeInfo created on the first traversal of the tree
 * for the corresponding parent
 * 
 * s: index into octree.nodes[] for the location of first non-NULL child
 * of parent
 *
 **********************************************************************/
void CreateSubtree(Node * parent, NodeInfo * parentInfo, int s){
    int i, childS_values[MaxChildren];    
    
    //setup table of s indexes for children
    for(i=0; i<MaxChildren; i++){
        childS_values[i] = s;
        if(parentInfo->childInfo[i]!=NULL){
            s+= parentInfo->childInfo[i]->subtree_size;
        }
    }
    
    //initialize and recurse
    for(i=0; i < MaxChildren; i++){
        if(parentInfo->childInfo[i] == NULL){
            parent->child[i] = NULL;             
        }
        else{
            #pragma omp task firstprivate(i)
            {
                parent->child[i] = InitializeNode(i, parentInfo->level + 1, parent, 
                    parentInfo->childInfo[i]->mid_x, parentInfo->childInfo[i]->mid_y, 
                    parentInfo->childInfo[i]->mid_z, parentInfo->childInfo[i]->low,
                    parentInfo->childInfo[i]->high, childS_values[i], 
                    parentInfo->childInfo[i]->isParent);   
                CreateSubtree(parent->child[i], parentInfo->childInfo[i],childS_values[i]+1);
            }
        }        
    }
    #pragma omp taskwait
}

/***********************************************************************
 * 
 * This is the topmost call for building tree data.
 * Create root_info and then call the main recursive funtion,
 * PartitionInfo(NodeInfo *)
 * 
 **********************************************************************/
NodeInfo * AnalyzeTree(int numParticles, int minimum_extent){
     
    nextInfo = 0;
    root_info = &info_array[nextInfo++];
    root_info->level = 0;
    root_info->low = 0;
    root_info->high = numParticles - 1;
    root_info->mid_x = minimum_extent+octree.edge_length[0]*0.5;
    root_info->mid_y = minimum_extent+octree.edge_length[0]*0.5;
    root_info->mid_z = minimum_extent+octree.edge_length[0]*0.5;
    root_info->isParent = False;
    PartitionInfo(root_info);    
    return root_info;    
}

/***********************************************************************
 * 
 * On a call to PartitionInfo(NodeInfo*), level, mid values and 
 * interval should already be set in the @info. The pointers to 
 * children info will be figured out and set by this function
 * using existing data in info
 * 
 **********************************************************************/
void PartitionInfo(NodeInfo * info){
    int c;
    double bisect;
    
    info->subtree_size = 1;
    for(c = 0; c < MaxChildren; c++){
        info->childInfo[c] = NULL;
    }
    //this sets precedent for defining a parent
    if((info->high - info->low + 1) > octree.maxBodiesPerNode) {
        if(maxDepth <= (info->level + 1)) //increment to lowest point of tree
        {
            maxDepth = info->level + 2;
        }
        info->isParent = True;
        //call subdivide
        subdivide(info->low, info->high, info->mid_x, info->mid_y, info->mid_z, info->childInts);
        bisect=octree.edge_length[info->level+1]*0.5;   
        if((bisect > 256)||(bisect<0)){
            (void)fprintf(stderr,"Bisect out of bounds: %f\n", bisect);
            exit(4);
        }
        //set children info pointers, recursing as needed
        for(c=0; c<MaxChildren; c++)
        {
            //-1 was a marker set in subdivide --> size of 0
            #pragma omp task firstprivate(c,info)
            {
                if(nextInfo >= floor(octree.nodeInfoCoeff*octree.nodeInfoDistribution))
                {
                    octree.threadFailures[threadId] = True;
                    //recursion killed by not calling PartitionInfo(...);

                    /*
                     * errouneously assign a pointer so parent doesn't break
                     * doesn't matter as tree will be rebuilt
                     */
                    info->childInfo[c] = &info_array[nextInfo-1];
                }
                else if(octree.threadFailures[threadId])
                {
                    //kills recursion of any tasks this thread picks up
                    //after the task that had original failure
                    /*
                     * errouneously assign a pointer so parent doesn't break
                     * doesn't matter as tree will be rebuilt
                     */
                    info->childInfo[c] = &info_array[nextInfo-1];
                }
                else
                {
                    NodeInfo * child_info = &info_array[nextInfo++];
                    child_info->low = info->childInts[c].low;
                    child_info->high = info->childInts[c].high;
                    //printf("Low: %i, High: %i\n", child_info->low, child_info->high);
                    child_info->level = info->level + 1;
                    child_info->mid_x = 
                        info->mid_x+(c & 1 ? bisect : -bisect);
                    child_info->mid_y = 
                        info->mid_y+(c & 2 ? bisect : -bisect);
                    child_info->mid_z = 
                        info->mid_z+(c & 4 ? bisect : -bisect);
                    child_info->isParent = False;
                    info->childInfo[c] = child_info;
                    PartitionInfo(child_info);
                }
            }
        }
        #pragma omp taskwait
        for(c = 0; c < MaxChildren; c++)
        {
            info->subtree_size+=info->childInfo[c]->subtree_size;
        }
     }
}

/***********************************************************************
 * 
 * Subdivides an existing octree, represented by its mid points
 * Fills in the passe in array of intervals corresponding to 8 children
 * 
 * This is really just a driver for the partitioning of bodies
 * 
 **********************************************************************/
void subdivide(int low, int high, float mid_x, float mid_y, float mid_z, Interval * childrenInts){
    int zDiv, yDiv1, yDiv2, xDiv1, xDiv2, xDiv3, xDiv4, i; 
    /*
        partition in this order:z then y then x. This order will guarantee that
        when child id's are tied to the subcubes, the ith cube represented
        by the partitioned array will correspond to thisParent->child[i] in the
        spatial sense. 
    */

    //z
    zDiv = PartitionBodies(low,high, mid_z, 2);
    //y
    if(zDiv<low)  //then none below
    {
        yDiv1 = low - 1; 
        yDiv2 = PartitionBodies(zDiv+1,high, mid_y, 1);
    }
    else if(zDiv == high) //none above
    {
        yDiv1 = PartitionBodies(low,zDiv, mid_y, 1);
        yDiv2 = high + 1;
    }
    else
    {
        yDiv1 = PartitionBodies(low,zDiv, mid_y, 1);
        yDiv2 = PartitionBodies(zDiv+1,high, mid_y, 1);
    }
    //x
    int intermed[5] = {low,yDiv1, zDiv, yDiv2, high};
    int *values[4] = {&xDiv1, &xDiv2, &xDiv3, &xDiv4};

    int currentLow = low;
    for(i=1; i<=4; i++)
    {
        if((intermed[i]>=currentLow)&&(intermed[i]<=high))
        {
            *values[i-1] = PartitionBodies(currentLow, intermed[i], mid_x, 0);
            currentLow = intermed[i]+1;
        }
        else{
            *values[i-1] = low -1;
        }
    }

    /*
        Now the partitions for sub-cubes 0 through 7 are:
        (low,xDiv1), (xDiv1,yDiv1), (yDiv1, xDiv2), (xDiv2, zDiv),
        (zDiv,xDiv3), (xDiv3,yDiv2), (yDiv2,xDiv4), (xDiv4, high)
    */
  
    int divisions[9] = {low, xDiv1, yDiv1, xDiv2, zDiv, xDiv3, yDiv2, xDiv4, high};
  
    for(i=0;i<MaxChildren;i++)
    {
        childrenInts[i].low = -1;
        childrenInts[i].high = -2;
    }
    int localLow = low;
    
    int COUNT = 0;
    for(i=0; i<MaxChildren; i++)
    {
        if((divisions[i+1]>=localLow)&&(divisions[i+1]<=high))
        {
            COUNT+= divisions[i+1]-localLow+1;
            childrenInts[i].low = localLow;
            childrenInts[i].high = divisions[i+1];
            localLow = divisions[i+1]+1;
            if(localLow>high){
                break;
            }
        }
    }
    //assert(COUNT==(high - low + 1));
}
void RebuildTree(unsigned long number_particles, int precision, double maximum_extent, double minimum_extent)
{
    double sTime;
    
    printf("************ Rebuilding tree: ");
    sTime = wcTime();
    
    FreeOctree();
    CreateOctree(number_particles, precision, maximum_extent, minimum_extent);
    
    printf(" %f\n", wcTime() - sTime);
    
}
int ReSort(Node *node, NodeInfo *info)
{
    int i, rebuild = False;
    if(node->isParent)
    {
        subdivide(node->pArrayLow, node->pArrayHigh, node->mid_x, node->mid_y, node->mid_z, info->childInts);
        for(i=0;i<MaxChildren;i++)
        {
            #pragma omp task firstprivate(i)
            {
                node->child[i]->pArrayLow = info->childInts[i].low;
                node->child[i]->pArrayHigh = info->childInts[i].high;
                if(ReSort(node->child[i], info->childInfo[i]))
                    rebuild = True;                    
            }
            
        }
        #pragma omp taskwait
    }
    else
    {
        rebuild = (node->pArrayHigh - node->pArrayLow + 1) > 2*octree.maxBodiesPerNode;
    }
    return rebuild;
}

/***********************************************************************
 * 
 * Setup thread private vars. Must be called directly within a 
 * parallel region
 * 
 **********************************************************************/
void threadSetup(){
    
    //increment and check function
    nextInfo = 0;
    maxDepth = 1;
    info_array = (NodeInfo*) malloc(floor(octree.nodeInfoCoeff*octree.nodeInfoDistribution)*sizeof(NodeInfo));
    if (info_array == (NodeInfo *) NULL) {
        fprintf(stderr, "unable to allocate memory for info_array\n");
        exit(4);
    }
    threadId = omp_get_thread_num();
    octree.threadFailures[threadId] = False;
}

/***********************************************************************
 * 
 * Setup lookup tables and other tree values
 * 
 **********************************************************************/
void treeVarSetup(int number_particles, int precision, double maximum_extent, double minimum_extent){
    double factorial[MaxPow+1], sum;
    float *ijk, *ijk_binomial;
    int i, j, k, n;
    
    //printf("Particles: %i\n",number_particles);
    octree.precision=precision;
    if (octree.precision > MaxPow)
        octree.precision=MaxPow;
    //printf("Precision: %lu\n",octree.precision);
    octree.coefficients=Binomial((int) octree.precision+3,3);
    //printf("Coefficients: %lu\n",octree.coefficients);
    /*
        ijk factorial
    */
    factorial[0]=1.0;
    sum=1.0;
    for (i=1; i <= octree.precision; i++)
    {
        sum*=(double) i;
        factorial[i]=sum;
    }
    octree.ijk_factorial=(float *) malloc(octree.coefficients*sizeof(float));
    if (octree.ijk_factorial == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }
    ijk=octree.ijk_factorial;
    for (i=0; i <= octree.precision; i++)
        for (j=0; j <= (octree.precision-i); j++)
            for (k=0; k <=  (octree.precision-i-j); k++)
                *ijk++=1.0/(factorial[i]*factorial[j]*factorial[k]);
    /*
        Ijk binomial used in constructing the A in formula 3.38
    */
    factorial[0]=1.0;
    sum=1.0;
    for (i=1; i <= octree.precision; i++)
    {
        sum*=(-0.5-(double) i+1.0);
        factorial[i]=sum;
    }
    octree.ijk_binomial=(float *) malloc(octree.coefficients*sizeof(float));
    if (octree.ijk_binomial == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }
    ijk_binomial=octree.ijk_binomial;
    ijk=octree.ijk_factorial;
    for (i=0; i <= octree.precision; i++)
        for (j=0; j <= (octree.precision-i); j++)
            for (k=0; k <= (octree.precision-i-j); k++)
                *ijk_binomial++=factorial[i+j+k]*(*ijk++);
    /*
        binomial
    */
    for (n=0; n <= octree.precision; n++)
        for (k=0; k <= octree.precision; k++)
            octree.binomial[n][k]=(float) Binomial(n,k);
    /*
        Array lookup
    */
    for (i=0; i <= (octree.precision+2); i++)
        octree.ref_three[i]=octree.coefficients-Binomial((int) octree.precision-i+2,3);
    for (i=0; i <= (octree.precision+2); i++)
        octree.ref_two[i]=Binomial((int) octree.precision-i+2,2);

    octree.depth = 1;

    octree.edge_length[0]=maximum_extent-minimum_extent;
    for (i=1; i < MaxTreeDepth; i++)
    {
        octree.edge_length[i]=octree.edge_length[i-1]*0.5;
    }
}

/***********************************************************************
 * 
 * Assigns important values to a node including coordicates, pointer to
 * parent, id, etc. 
 * 
 **********************************************************************/
Node *InitializeNode(int id,int level,Node *parent, double mid_x, 
    double mid_y,double mid_z, int low, int high, int index, 
    int isParent)
{
  int i;
  
  Node * nodeToInit =  &octree.nodes[index];
  
  nodeToInit->parent=parent;
  
  for (i=0; i < MaxChildren; i++)
    nodeToInit->child[i]=(Node *) NULL;
  nodeToInit->id=id;
  nodeToInit->level=level;
  nodeToInit->mid_x=mid_x;
  nodeToInit->mid_y=mid_y;
  nodeToInit->mid_z=mid_z;
  nodeToInit->isParent = isParent;
  nodeToInit->pArrayLow = low;
  nodeToInit->pArrayHigh = high;
  for(i=0;i<4;i++)
  {
    nodeToInit->phi[i]= &(octree.PHIs[i][octree.coefficients*index]);
    nodeToInit->psi[i]= &(octree.PSIs[i][octree.coefficients*index]);
  }
  
  //marker to tell if di node or not
  nodeToInit->DI_tag = -1;

    //inefficient
  for(i=0; i<MaxList1; i++){
    nodeToInit->Colleagues[i]=NULL;
    nodeToInit->list1[i]=NULL;
  }
  for(i = 0; i < MaxList3; i++){
    nodeToInit->list2[i]=NULL;
    nodeToInit->list3[i]=NULL;
    nodeToInit->list4[i]=NULL;
  }
  
    nodeToInit->colleagueCount = &octree.counts[5*index];    
    nodeToInit->list1Count = &octree.counts[5*index+1];
    nodeToInit->list2Count = &octree.counts[5*index+2];  
    nodeToInit->list4Count = &octree.counts[5*index+3];  
    nodeToInit->list3Count = &octree.counts[5*index+4];

  return nodeToInit;
}

/***********************************************************************
*
*  ParticleSwap - a deep copy that utilizes the fact that we want to 
*  swap the contents of the two particles handed in. Swap in
*  octree.bodies
* 
***********************************************************************/
void ParticleSwap(int from, int to)
{
  Particle temp;
  temp.position = octree.bodies[from].position;
  temp.force = octree.bodies[from].force;
  
  //copy data
//   temp.position[0] = octree.bodies[from].position[0];
//   temp.position[1] = octree.bodies[from].position[1];
//   temp.position[2] = octree.bodies[from].position[2];
//   
//   temp.force[0] = octree.bodies[from].force[0];
//   temp.force[1] = octree.bodies[from].force[1];
//   temp.force[2] = octree.bodies[from].force[2];
//   
// 
//   //overwrite
  octree.bodies[from].position = octree.bodies[to].position;
  octree.bodies[from].force = octree.bodies[to].force;
//   octree.bodies[from].position[0] = octree.bodies[to].position[0];
//   octree.bodies[from].position[1] = octree.bodies[to].position[1];
//   octree.bodies[from].position[2] = octree.bodies[to].position[2];
//   
//   octree.bodies[from].force[0] = octree.bodies[to].force[0];
//   octree.bodies[from].force[1] = octree.bodies[to].force[1];
//   octree.bodies[from].force[2] = octree.bodies[to].force[2];
//   
// 
//   //copy over
  octree.bodies[to].position = temp.position;
  octree.bodies[to].force = temp.force;
//   octree.bodies[to].position[0] = temp.position[0];
//   octree.bodies[to].position[1] = temp.position[1];
//   octree.bodies[to].position[2] = temp.position[2];
//   
//   octree.bodies[to].force[0] = temp.force[0];
//   octree.bodies[to].force[1] = temp.force[1];
//   octree.bodies[to].force[2] = temp.force[2];

}

/***********************************************************************
 * 
 * PartitionBodies() partitions the particles array up by the x, y, 
 * and z coordinates of each particle based on passed in mid value. Only 
 * does one dimension of partitioning per called, which dimension is 
 * determined by 'axis' parameter.
 * 
 **********************************************************************/
int PartitionBodies(int low, int high, float midValue, int axis)
{
    int zeroTop, j;
  
    zeroTop = low -1; 
    //assert((axis==0) || (axis==1) || (axis==2));

    for (j = low; j<= high; j++)
    {
      if(octree.bodies[j].position[axis] <= midValue)
      {
        ParticleSwap(++zeroTop, j);
      }
    }

    return zeroTop;
}

void constructTreeInfo(int numParticles, int minimum_extent){
    #pragma omp parallel
    {
        #pragma omp critical
        {
            threadSetup();
        }
        //top level call to build tree of node info
        #pragma omp single
        {
            octree.rootInfo = AnalyzeTree(numParticles, minimum_extent);
        }        
    }
}

/***********************************************************************
 *
 * ParPrefixSum()
 * @brief - in place
 *
 **********************************************************************/
void ParPrefixSum(unsigned int *A, unsigned int n)
{
    unsigned int i, j, f;

    #pragma omp parallel
    {
        octree.numThreads = omp_get_num_threads();
    }

    unsigned int b[octree.numThreads + 1], PrevSums[octree.numThreads];

    f = (unsigned int) floor(n/octree.numThreads);
    for(i=0; i<=octree.numThreads; i++)
    {
        b[i] = i*f;
    }
    b[octree.numThreads] = n;

    #pragma omp parallel
    {
        #pragma omp single
        {
            for(i=0; i<octree.numThreads; i++)
            {
                #pragma omp task firstprivate(i,j)
                {
                    for(j=b[i]+1; j<b[i+1]; j++)
                    {
                        A[j] += A[j - 1];
                    }
                }
            }
        }
        #pragma omp barrier
    }

    PrevSums[0] = 0;
    for(i=1;i<octree.numThreads;i++)
    {
        PrevSums[i] = PrevSums[i-1] + A[b[i]-1];
    }

    #pragma omp parallel
    {
        #pragma omp single
        {
            for(i=1; i<octree.numThreads; i++)
            {
                #pragma omp task firstprivate(i,j)
                {
                    for(j=b[i]; j<b[i+1]; j++)
                    {
                        A[j] += PrevSums[i];
                    }
                }
            }
        }
        #pragma omp barrier
    }
}
