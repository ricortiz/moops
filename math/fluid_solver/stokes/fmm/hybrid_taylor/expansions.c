#include "hybridfmm.h"
#include <assert.h>
/***********************************************************************
 * NOTE 1:
 *
 * A simplification is used throughout this code. Note that if a outer
 * to inner shift is performed as described by Lemma 3.2.2 on pg 23 then
 * each "c" coefficient has a 1/i!j!k! on the front of it as defined by
 * formula 3.35. Now if we shift this expansion of the form of formula
 * 3.34 down to a child as defined by Lemma 3.2.3 then as defined in
 * formula 3.42, each of these "c" coefficients will be multiplied by
 * the binomial coefficients in formula 3.42 ( (i+a) choose (a),
 * (j+b) choose (b) and (k+g) choose (g)). Note that (1/(i!j!k!)) * the
 * binomial coefficients just listed = (1/(i!j!k!)) * (1/(a!b!g!)).
 *
 * So if we know that we are going to downshift an expansion we form
 * from Lemma 3.2.2 of the form 3.34, then if we leave off the 1/i!j!k!
 * on each of the "c" coefficients in 3.35 then when we perform the
 * downshift and carry out formula 3.42 we do not need to calculate the
 * binomial coefficients in 3.42. Instead we substitue these for 1/a!b!g!
 * and then in addition multiply by 1/i!j!k!.
 *
 * Now we can only do then when we know we will later perform the
 * downshift. This is true for all interior nodes. Hence for all
 * interior nodes, we can skip multiplying by 1/i!j!k! in 3.35. And then
 * when we evaluate the downshift in 3.42 we do not have to evaluate the
 * binomial coefficients. However for leaf nodes we must multiply by
 * 1/i!j!k! in 3.35 since we do not perform downshifts on lead nodes.
 *
 * This change describes the variations throughout the expansion
 * functions from the methods defined in the paper. The advantage to
 * this is that we can reuse octree.ijk_factorial which is the 1/i!j!k!
 * values when we evaluate 1/a!b!g! as opposed to having to create
 * additional lookup tables for the binomial coefficients of formula 3.42
 *
 * To understand this note that when we shift into the
 * inner expansion of an interior node ShiftOuterToInnerForNonLeaf()
 * we did not multiply by 1/i!j!k! as
 *
 *
 **********************************************************************/

/***********************************************************************
 *
 * FormOuterExpansion(Node *node)
 * @brief - used in the upsweep phase. Given a node computes the
 * multipole expansion valid at distance outside the node which
 * represents the potential field due to the particles in that node.
 * The phi coefficients computed are the "a" coefficients of formula 3.7
 * on pg 14
 *
 **********************************************************************/
void FormOuterExpansion(Node *node)
{
    assert(node->pArrayLow >= 0);
    double x_distance, y_distance, z_distance, sum;
    unsigned long i, j, k, n;

    Particle *p;
    int low, high, pNum, l;
    float x_power[MaxPow+1], y_power[MaxPow+1], z_power[MaxPow+1], *ijk,
    *phi;

    low = node->pArrayLow;
    high = node->pArrayHigh;

    /*
     *     Form multipole expansions of potential field due to particles
     *     in each cube about the cube center at the leaf nodes.
     */
    for (l = 0;l < 4;l++)

    {
        for (pNum = low; pNum <= high; pNum++)
        {
            p = &octree.bodies[pNum];
            x_distance = p->position[0] - node->mid_x + Epsilon;
            y_distance = p->position[1] - node->mid_y + Epsilon;
            z_distance = p->position[2] - node->mid_z + Epsilon;
            for (i = 0; i <= octree.precision; i++)
            {
                x_power[i] = pow(x_distance, (double) i);
                y_power[i] = pow(y_distance, (double) i);
                z_power[i] = pow(z_distance, (double) i);
            }

            float fdotx = p->position[0] * p->force[0] + p->position[1] * p->force[1] + p->position[2] * p->force[2];
            float charges[4] = {p->force[0], p->force[1], p->force[2], fdotx};
            phi = node->phi[l];
            ijk = octree.ijk_factorial;
            for (i = 0; i <= octree.precision; i++)
            {
                for (j = 0; j <= (octree.precision - i); j++)
                {
                    for (k = 0; k <= (octree.precision - i - j); k++)
                    {
                        *phi += charges[l] * x_power[i] * y_power[j] * z_power[k] * (i + j + k & 1 ? -1 : 1) * (*ijk++);
                        phi++;
                    }
                }
            }
        }
    }
}

/***********************************************************************
 *
 * ShiftFromChildToParent(Node *node)
 * @brief - to be called only on parent nodes. For each child of @node
 * the outer expansions of each child is combined into one expansion
 * centered at the center of @node. This expansion represents
 * contributions from all particles in @node
 *
 * This performs the conversion specified by Lemma 3.2.1 pg 21. The "a"
 * coefficients in formula 3.30 are the phi coefficients for the child
 * the "a'" coefficients are formed by calling MultipoleExpansion
 *
 **********************************************************************/
void ShiftFromChildToParent(Node *node)
{
    assert(node->pArrayLow >= 0);
    double sum;
    unsigned long i, j, k, a, b, g;
    float *phi;
    int id, l;

    //temp coefficients
    float phi_tilde[octree.coefficients];

    for (l = 0;l < 4;++l)
    {
        /*
         *           for each child of this node, which is a parent by if statement,
         *           shift the child's phi to parent and add to parent
         */
        for (id = 0; id < 8; id++)
        {
            if (node->child[id]->pArrayLow != -1)
            {
                MultipoleExpansion(node->child[id], node, phi_tilde);
                phi = node->phi[l];
                for (i = 0; i <= (long) octree.precision; i++)
                    for (j = 0; j <= (long) (octree.precision - i); j++)
                        for (k = 0; k <= (long) (octree.precision - i - j); k++)
                        {
                            for (a = 0; a <= i; a++)
                                for (b = 0; b <= j; b++)
                                    for (g = 0; g <= k; g++)
                                    {
                                        *phi += ArrayLookup(node->child[id]->phi[l], a, b, g) *
                                                ArrayLookup(phi_tilde, i - a, j - b, k - g);
                                    }
                            phi++;
                        }
            }
        }
    }
}


/***********************************************************************
 *
 * OuterToInner()
 * @brief - this function shifts an outer expansion which represents
 * the potential field due to the particles in the @source node to an
 * inner expansion centered at the @target node. This is function
 * performs Lemma 3.2.2 found on page 23.
 *
 * Depending on whether or not the target is a leaf node, multiply by
 * 1/i!j!k!. This differs from the definition in the paper, for an
 * explaination of this see Note 1 at top of page.
 *
 **********************************************************************/
void OuterToInner(Node *source, Node *target, float *target_psi[4], int isLeaf)
{
    assert((source->pArrayLow >= 0) && (target->pArrayLow >= 0));
    int i, j, k, n, a, b, g, l;
    float psi_tilde[octree.coefficients], *psi, *phi, *ijk;
    double sum;

    for (l = 0;l < 4;++l)
    {
        /*
         * Dump some coefficients into psi_tilde
         * See FormLocalExpansion for description of coefficients
         */
        FormLocalExpansion(source, target, psi_tilde);
        psi = target_psi[l];
        phi = source->phi[l];
        ijk = octree.ijk_factorial;
        for (i = 0; i <= (long) octree.precision; i++)
            for (j = 0; j <= (long) (octree.precision - i); j++)
                for (k = 0; k <= (long) (octree.precision - i - j); k++)
                {
                    sum = 0.0;
                    n = i + j + k;
                    for (a = 0; a <= (long) (octree.precision - n); a++)
                    {
                        for (b = 0; b <= (long) (octree.precision - n - a); b++)
                        {
                            for (g = 0; g <= (long) (octree.precision - n - a - b); g++)
                            {
                                sum += ArrayLookup(phi, a, b, g) *
                                       ArrayLookup(psi_tilde, i + a, j + b, k + g);

                            }
                        }
                    }
                    if (isLeaf)
                        *psi += sum * (*ijk++);
                    else
                        *psi += sum;
                    psi++;
                }
    }
}

/***********************************************************************
 *
 * DownShift()
 * @brief - this function performs a "downshift" as defined by Lemma
 * 3.2.3 page 25. The "c" coefficients are the psi coefficients of the
 * @Parent. Note that this is not in the exact form as the paper defines
 * in 3.42 as no binomial coefficients are calculated. For explaination
 * of this see Note 1 at the top of the page.
 *
 * Depending on whether or not the target is a leaf node, multiply by
 * 1/i!j!k!. This differs from the definition in the paper, for an
 * explaination of this see Note 1 at top of page.
 *
 **********************************************************************/
void DownShift(Node *Child, int isLeaf)
{
    assert(Child->pArrayLow >= 0);
    float x_power[MaxPow+1], y_power[MaxPow+1], z_power[MaxPow+1], *ijk,

    *psi;
    double x_distance, y_distance, z_distance, sum;
    int i, j, k, n, a, b, g, l;
    Node *Parent = Child->parent;

    x_distance = Child->mid_x - Parent->mid_x + Epsilon;
    y_distance = Child->mid_y - Parent->mid_y + Epsilon;
    z_distance = Child->mid_z - Parent->mid_z + Epsilon;
    for (i = 0; i <= (long) octree.precision; i++)

    {
        x_power[i] = pow(x_distance, (double) i);
        y_power[i] = pow(y_distance, (double) i);
        z_power[i] = pow(z_distance, (double) i);
    }

    for (l = 0;l < 4;++l)
    {
        psi = Child->psi[l];
        ijk = octree.ijk_factorial;
        for (i = 0; i <= (long) octree.precision; i++)
            for (j = 0; j <= (long) (octree.precision - i); j++)
                for (k = 0; k <= (long) (octree.precision - i - j); k++)
                {
                    sum = 0.0;
                    n = i + j + k;
                    for (a = 0; a <= (long) (octree.precision - n); a++)
                        for (b = 0; b <= (long) (octree.precision - n - a); b++)
                            for (g = 0; g <= (long) (octree.precision - n - a - b); g++)
                            {
                                sum += ArrayLookup(Parent->psi[l], i + a, j + b, k + g) *
                                       ArrayLookup(octree.ijk_factorial, a, b, g) *
                                       x_power[a] * y_power[b] * z_power[g];
                            }
                    if (isLeaf)
                        *psi += sum * (*ijk++);
                    else
                        *psi += sum;
                    psi++;
                }
    }
}

/***********************************************************************
 *
 * ApplyLocalExpansion(Node *node)
 * @brief - to be called on leaf nodes. Evaluates @node's local
 * expansion at every particle in the node
 *
 **********************************************************************/
void ApplyLocalExpansion(Node *node)
{
    assert(node->pArrayLow >= 0);
    double dx, dy, dz;
    Particle *p;
    int i, j, k, pNum;
    float x_power[MaxPow+1], y_power[MaxPow+1], z_power[MaxPow+1], *psi_0,
    *psi_1, *psi_2, *psi_3, dx_power[MaxPow+1], dy_power[MaxPow+1],

    dz_power[MaxPow+1];

    for (pNum = node->pArrayLow; pNum <= node->pArrayHigh; pNum++)
    {
        p = &octree.bodies[pNum];

        dx=p->position[0]-node->mid_x;
        dy=p->position[1]-node->mid_y;
        dz=p->position[2]-node->mid_z;

        x_power[0]=1;
        y_power[0]=1;
        z_power[0]=1;

        dx_power[0] = 0;
        dy_power[0] = 0;
        dz_power[0] = 0;
        for (i=1; i <= octree.precision; i++)
        {
            x_power[i]=pow(dx,(double) i);
            y_power[i]=pow(dy,(double) i);
            z_power[i]=pow(dz,(double) i);

            dx_power[i] = i*pow(dx,(double) i-1.0);
            dy_power[i] = i*pow(dy,(double) i-1.0);
            dz_power[i] = i*pow(dz,(double) i-1.0);

            if (isnanf(dx_power[i])||isnanf(dy_power[i])||isnanf(dz_power[i])
                    ||isnanf(x_power[i])||isnanf(y_power[i])||isnanf(z_power[i]))
            {
                printf("****************NaN\n");
                exit(4);
            }
        }
        psi_0 = node->psi[0];
        psi_1 = node->psi[1];
        psi_2 = node->psi[2];
        psi_3 = node->psi[3];
        for (i = 0; i <= octree.precision; i++)
        {
            for (j = 0; j <= (octree.precision - i); j++)
            {
                for (k = 0; k <= (octree.precision - i - j); k++)

                {
                    octree.potentials[pNum].x += (*psi_0) * x_power[i] * y_power[j] * z_power[k];
                    octree.potentials[pNum].y += (*psi_1) * x_power[i] * y_power[j] * z_power[k];
                    octree.potentials[pNum].z += (*psi_2) * x_power[i] * y_power[j] * z_power[k];

                    octree.fields[pNum].field[0][0] += (*psi_0) * dx_power[i] * y_power[j] * z_power[k];
                    octree.fields[pNum].field[1][0] += (*psi_0) * x_power[i] * dy_power[j] * z_power[k];
                    octree.fields[pNum].field[2][0] += (*psi_0) * x_power[i] * y_power[j] * dz_power[k];

                    octree.fields[pNum].field[0][1] += (*psi_1) * dx_power[i] * y_power[j] * z_power[k];
                    octree.fields[pNum].field[1][1] += (*psi_1) * x_power[i] * dy_power[j] * z_power[k];
                    octree.fields[pNum].field[2][1] += (*psi_1) * x_power[i] * y_power[j] * dz_power[k];

                    octree.fields[pNum].field[0][2] += (*psi_2) * dx_power[i] * y_power[j] * z_power[k];
                    octree.fields[pNum].field[1][2] += (*psi_2) * x_power[i] * dy_power[j] * z_power[k];
                    octree.fields[pNum].field[2][2] += (*psi_2) * x_power[i] * y_power[j] * dz_power[k];

                    octree.fields[pNum].field[0][3] += (*psi_3) * dx_power[i] * y_power[j] * z_power[k];
                    octree.fields[pNum].field[1][3] += (*psi_3) * x_power[i] * dy_power[j] * z_power[k];
                    octree.fields[pNum].field[2][3] += (*psi_3) * x_power[i] * y_power[j] * dz_power[k];

                    psi_0++;
                    psi_1++;
                    psi_2++;
                    psi_3++;
                }
            }
        }
        unsigned int idx = 3 * pNum;
        octree.CPU_Veloc[idx]   +=  0.039788735772974 * (octree.potentials[pNum].x - p->position[0] * octree.fields[pNum].field[0][0] - p->position[1] * octree.fields[pNum].field[0][1] - p->position[2] * octree.fields[pNum].field[0][2] + octree.fields[pNum].field[0][3]);
        octree.CPU_Veloc[idx+1] +=  0.039788735772974 * (octree.potentials[pNum].y - p->position[0] * octree.fields[pNum].field[1][0] - p->position[1] * octree.fields[pNum].field[1][1] - p->position[2] * octree.fields[pNum].field[1][2] + octree.fields[pNum].field[1][3]);
        octree.CPU_Veloc[idx+2] += 0.039788735772974 * (octree.potentials[pNum].z - p->position[0] * octree.fields[pNum].field[2][0] - p->position[1] * octree.fields[pNum].field[2][1] - p->position[2] * octree.fields[pNum].field[2][2] + octree.fields[pNum].field[2][3]);

    }
}

/***********************************************************************
 *
 * FormLocalExpansion
 * @brief - this function is used for shifting an outer expansion to
 * an inner expansion. This function fills @psi_tilde with the product
 * of the "b" coefficients and (i+a)!, (j+b)! and (k + g)! found in
 * formula 3.35, pg 24.
 *
 * @note : this implementation does not formulate the b's exactly as
 * 3.38 specifies. See Note 1 at top of page
 *
 **********************************************************************/
void FormLocalExpansion(Node *interactive_neighbor, Node *node, float *psi_tilde)
{
    assert((interactive_neighbor->pArrayLow >= 0) && (node->pArrayLow >= 0));
    double distance, tau_x, tau_y, tau_z, dx, dy, dz, sum;
    long i, j, k, n, a, b, g;
    float r_zero[MaxPow+1], x_power[MaxPow+1], y_power[MaxPow+1],

    z_power[MaxPow+1], *ijk;

    dx = interactive_neighbor->mid_x - node->mid_x + Epsilon;
    dy = interactive_neighbor->mid_y - node->mid_y + Epsilon;
    dz = interactive_neighbor->mid_z - node->mid_z + Epsilon;
    distance = sqrt(dx * dx + dy * dy + dz * dz);
    tau_x = 2.0 * (dx / distance);
    tau_y = 2.0 * (dy / distance);
    tau_z = 2.0 * (dz / distance);
    for (i = 0; i <= (long) octree.precision; i++)
    {
        x_power[i] = pow(tau_x, (double) i);
        y_power[i] = pow(tau_y, (double) i);
        z_power[i] = pow(tau_z, (double) i);
        r_zero[i] = 1.0 / pow(distance, (double) i + 1);
    }
    ijk = octree.ijk_factorial;
    for (i = 0; i <= (long) octree.precision; i++)
    {
        for (j = 0; j <= (long) (octree.precision - i); j++)
        {
            for (k = 0; k <= (long) (octree.precision - i - j); k++)
            {
                sum = 0.0;
                for (a = 0; a <= (i / 2); a++)
                    for (b = 0; b <= (j / 2); b++)
                        for (g = 0; g <= (k / 2); g++)
                            sum += ArrayLookup(octree.ijk_binomial, i - a, j - b, k - g) *
                                   octree.binomial[i-a][a] * octree.binomial[j-b][b] *
                                   octree.binomial[k-g][g] * x_power[i-2*a] * y_power[j-2*b] *
                                   z_power[k-2*g];
                n = i + j + k;
                sum *= (n & 1 ? -1 : 1) * r_zero[n] / (*ijk++);
                *psi_tilde++ = sum;

            }
        }
    }
}
/***********************************************************************
 *
 * MultipoleExpansion()
 * @brief - this function is used as part of shifting from a child to
 * its parent. This function fills in the "a'" coefficients in formula
 * 3.30 on pg 21. The "a'" coefficients are places into phi_tilde
 *
 **********************************************************************/
void MultipoleExpansion(Node *node, Node *parent, float *phi_tilde)
{
    assert((parent->pArrayLow >= 0) && (node->pArrayLow >= 0));

    double dx, dy, dz;
    long i, j, k;
    float x_power[MaxPow+1], y_power[MaxPow+1], z_power[MaxPow+1], *ijk;

    dx = node->mid_x - parent->mid_x + Epsilon;
    dy = node->mid_y - parent->mid_y + Epsilon;
    dz = node->mid_z - parent->mid_z + Epsilon;
    for (i = 0; i <= (long) octree.precision; i++)
    {
        x_power[i] = pow(-dx, (double) i);
        y_power[i] = pow(-dy, (double) i);
        z_power[i] = pow(-dz, (double) i);
    }
    ijk = octree.ijk_factorial;
    for (i = 0; i <= (long) octree.precision; i++)
        for (j = 0; j <= (long) (octree.precision - i); j++)
            for (k = 0; k <= (long) (octree.precision - i - j); k++)
                *phi_tilde++ = x_power[i] * y_power[j] * z_power[k] * (*ijk++);
}


