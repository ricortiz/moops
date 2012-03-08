
extern "C" {
    #include "hybridfmm.h"
}
/***********************************************************************
 *
 * Main
 *
 **********************************************************************/
int main(int argc, char **argv)
{
    double maximum_extent, minimum_extent, startTime, secondaryTime,
           dt, tStep, fmmTime, cTime;
    int i, j, aflag = 0, pflag = 0, cflag = 0, dflag = 0, c, pCount = 0;
    unsigned long number_particles, precision;
    char *params[argc];
    unsigned int iterations;

    Particle *tmpBod;
    float *temp_GPU_velocities;

    for(i = 0; i < argc; i++)
    {
        if(strchr(argv[i], '-') == NULL)
            params[pCount++] = argv[i];
    }

    //fprintf(stderr, "--------------------------START RUN-------------------------\n");

    number_particles = 1000;
    if (pCount > 1)
    {
        number_particles = (unsigned long) atol(params[1]);
        octree.numParticles = (unsigned int) atoi(params[1]);
    }
    printf("Number of particles: %lu\n", number_particles);
    precision = 6;
    if (pCount > 2)
    {
        precision = (unsigned long) atol(params[2]);
    }
    printf("Precision: %lu\n", precision);
    octree.maxBodiesPerNode = 50;
    if (pCount > 3)
    {
        octree.maxBodiesPerNode = (unsigned long) atol(params[3]);
        if(octree.maxBodiesPerNode < 1)
        {
            printf("max bodies parameter must be greater than zero\n");
            exit(4);
        }
    }
    printf("S: %i\n", octree.maxBodiesPerNode);

    iterations = 1.0;
    if(pCount > 4)
    {
        iterations = (unsigned int) atoi(params[4]);
    }
    printf("Iterations: %i\n", iterations);

//     tmpBod = (Particle *) malloc((number_particles * sizeof(Particle)) + 2 * MEMORY_ALIGNMENT);
    tmpBod = new Particle[number_particles+ 1 << 11];
    if (tmpBod == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }

//     temp_GPU_velocities = (float *) malloc(number_particles * sizeof(float) + 2 * MEMORY_ALIGNMENT);
temp_GPU_velocities = new float[3*number_particles + 1 << 11];
    if (temp_GPU_velocities == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }
    //align for GPU
    octree.GPU_Veloc = (float*) ALIGN_UP( temp_GPU_velocities, MEMORY_ALIGNMENT );
    octree.bodies = (Particle *) ALIGN_UP( tmpBod, MEMORY_ALIGNMENT );


//     octree.CPU_Veloc = (double*) malloc(number_particles * sizeof(double/*CPU_Velocities*/));
    octree.CPU_Veloc = new float[3*number_particles];
    if (octree.CPU_Veloc == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }
    octree.fields = (Field *) malloc(number_particles * sizeof(Field));
    if (octree.fields == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }
    octree.potentials = (Potential *) malloc(number_particles * sizeof(Potential));
    if (octree.potentials == NULL)
    {
        printf("Not able to allocate memory");
        exit(4);
    }
    //---get options---//
    while ((c = getopt (argc, argv, "apcd")) != -1)
    {
        switch (c)
        {
            case 'a':
                aflag = 1;
                break;
            case 'p':
                pflag = 1;
                break;
            case 'c':
                cflag = 1;
                break;
            case 'd':
                dflag = 1;
                break;
            default:
                break;
        }
    }

    srand(30428877);
    minimum_extent = MinRange;
    maximum_extent = MaxRange;

    memset(octree.GPU_Veloc, 0.0, number_particles * sizeof(GPU_Velocities));
    memset(octree.CPU_Veloc, 0.0, number_particles * sizeof(CPU_Velocities));
    memset(octree.potentials, 0.0, number_particles * sizeof(Potential));
    memset(octree.fields, 0.0, number_particles * sizeof(Field));

    if (!pflag)
    {
        for (i = 0; i < (long) number_particles; i++)
        {
            octree.bodies[i].position[0] = (rand() * ((MaxRange - MinRange) + MinRange)) / RAND_MAX;
            octree.bodies[i].position[1] = (rand() * ((MaxRange - MinRange) + MinRange)) / RAND_MAX;
            octree.bodies[i].position[2] = (rand() * ((MaxRange - MinRange) + MinRange)) / RAND_MAX;

            octree.bodies[i].force[0] = .001;
            octree.bodies[i].force[1] = 0.0;
            octree.bodies[i].force[2] = 0.0;
        }
    }
    else
    {
        Plummer(number_particles);
    }

    octree.output = fopen("output.csv", "a");

    //create octree
    startTime = wcTime();
    CreateOctree(number_particles, precision, maximum_extent, minimum_extent);
    cTime = wcTime() - startTime;

    //move this into like variable setup
    /*octree.scale = (double *) malloc(octree.depth*sizeof(double));
    octree.scale[0] = 1.0/MaxRange;
    for (i = 1; i < octree.depth; i++ ) {
            octree.scale[i]  = 2.0*octree.scale[i-1];
        }*/

    //time step vars
    tStep = 0.001;
    dt = tStep / (double)iterations;

    printf("dt: %.15f\n", dt);

    //call allpairs or FMM
    if(aflag)
    {
        startTime = wcTime();
        for(i = 1; i <= iterations; i++)
        {
            AllPairs(number_particles, dt, .05);
        }
        printf("Allpairs: %8.3fs\n", wcTime() - startTime);
        //RebuildTree(number_particles, precision, maximum_extent, minimum_extent);
    }
    else
    {
        fmmTime = wcTime();
        for(i = 1; i <= iterations; i++)
        {
            fprintf(octree.output, "%lu, %d, %lu, ", number_particles, octree.maxBodiesPerNode, precision);
            FMM(dt);
            //reset phi and psi
            for(j = 0; j < 4; j++)
            {
                memset(octree.PHIs[j], 0.0, octree.rootInfo->subtree_size * octree.coefficients * sizeof(float));
                memset(octree.PSIs[j], 0.0, octree.rootInfo->subtree_size * octree.coefficients * sizeof(float));
            }
            memset(octree.potentials, 0.0, number_particles * sizeof(Potential));
            memset(octree.fields, 0.0, number_particles * sizeof(Field));
            //secondaryTime = wcTime();
            //if(!dflag)
            //    printf("Performing iteration %i ...  %.3fs\n",i, wcTime() - secondaryTime);
            //fprintf(octree.output, " %f\n", wcTime()-secondaryTime);
        }
        fmmTime = wcTime() - fmmTime;
        fprintf(octree.output, ", %f\n", fmmTime + cTime);
    }
    /*
    for(i=0;i<number_particles;i++)
    {
          fprintf(octree.output, "%.10f, %.10f, %.10f\n", octree.bodies[i].position[0], octree.bodies[i].position[1], octree.bodies[i].position[2]);
        }

        for(i=0;i<10;i++)
        {
              fprintf(stderr, "Forces x: %.10f, y: %.10f, z: %.10f\n", octree.CPU_Forces[i].x, octree.CPU_Forces[i].y, octree.CPU_Forces[i].z);
            }*/

    //Free up memory
    delete [] tmpBod;
    delete [] temp_GPU_velocities;
    delete [] octree.CPU_Veloc;
//     free((Field *) octree.fields);
//     free((double *) octree.potentials);
    FreeOctree();

    fclose(octree.output);
    //fprintf(stderr, "--------------------------END RUN-------------------------\n");
    return(False);
}