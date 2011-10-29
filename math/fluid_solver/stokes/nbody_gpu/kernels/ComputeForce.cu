#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "cuda.h"
#include "mex.h"

#include "Forces_Kernel.cu"

#define NUM_POINTS_TOWER 120

//prototypes
void PrintFromMatlab(FILE *, float *, float , float , float , float ,	int );
void Print2dArray(FILE *, char *,  float * , int , int );
void PrintToMatlab(FILE *, float *);
void PrintResults_StrengthKernel(FILE *fp, int * , float * , float * , float * , float *  , float *  , float *   , float *  , float *  );

//global 
FILE *fmex, *fcuda;
int ROWS, COLS;
int num_sources, num_towers;
int col_array[3] = {0,1,2};

//********************************************************* Predictor  Loop  ******************************************************************

extern "C"
{
	
	void 	ForceEval(float * h_sources,float *  h_strengths, float resting_length, float vertical_length, float diagonal_length, float spring_constant)
	{						
		
		//fcuda = fopen("logcuda.txt","a");
		
		//for logging - begin
		/*int *hlog_indices, *dlog_indices;
		float *hlog_positions, *hlog_strengths, *hlog_c, *hlog_pos1, *hlog_pos2, *hlog_magnitude, *hlog_springk, *hlog_length;
		float *dlog_positions, *dlog_strengths, *dlog_c,  *dlog_pos1, *dlog_pos2, *dlog_magnitude, *dlog_springk, *dlog_length;
		
		hlog_indices = new int [10];
		hlog_positions = new float [30];
		hlog_strengths = new float [30];
		hlog_c = new float [10];
		hlog_pos1 = new float [30];
		hlog_pos2 = new float [30];
		hlog_magnitude = new float [10];
		hlog_springk = new float [10];
		hlog_length = new float [10];
		
		memset(hlog_positions, 0, 30 * sizeof(float));
		memset(hlog_strengths, 0, 30 * sizeof(float));
		memset(hlog_c, 0, 10 * sizeof(float));
		memset(hlog_pos1, 0, 30 * sizeof(float));
		memset(hlog_pos2, 0, 30 * sizeof(float));
		memset(hlog_magnitude, 0, 10 * sizeof(float));
		memset(hlog_springk, 0, 10 * sizeof(float));
		memset(hlog_length, 0, 10 * sizeof(float));	
		
		cudaMalloc((void**)&dlog_indices, 10 * sizeof(int));   
		cudaMalloc((void**)&dlog_positions, 30 * sizeof(float));    
		cudaMalloc((void**)&dlog_strengths, 30 * sizeof(float)); 
		cudaMalloc((void**)&dlog_c, 10 * sizeof(float)); 
		cudaMalloc((void**)&dlog_pos1, 30 * sizeof(float)); 
		cudaMalloc((void**)&dlog_pos2, 30 * sizeof(float)); 
		cudaMalloc((void**)&dlog_magnitude, 10 * sizeof(float)); 
		cudaMalloc((void**)&dlog_springk, 10 * sizeof(float)); 
		cudaMalloc((void**)&dlog_length, 10 * sizeof(float)); */
		//for logging - end		
				
		int NUM_THREADS = NUM_POINTS_TOWER;
		int SIZE = num_sources*3*sizeof(float);
		float *d_sources, *d_strengths;
		
		cudaMalloc((void**)&d_sources, SIZE);   
		cudaMalloc((void**)&d_strengths, SIZE);    
		cudaMemcpy(d_strengths, h_strengths, SIZE, cudaMemcpyHostToDevice);
		cudaMemcpy(d_sources, h_sources, SIZE, cudaMemcpyHostToDevice);
			
		int blockSize = NUM_THREADS;
		int nBlocks = (int)num_sources/blockSize + (num_sources%blockSize == 0?0:1);	  
		Kernel_Eval_Strengths<<< nBlocks, blockSize >>> ((float3*)d_strengths,(float3*)d_sources, spring_constant ,	resting_length, vertical_length, diagonal_length, num_sources, num_towers);
																									//dlog_indices, dlog_strengths, dlog_positions , dlog_c, dlog_pos1, dlog_pos2, dlog_magnitude, dlog_springk, dlog_length);

		// copy memory from device to host
		cudaMemcpy(h_strengths, d_strengths, SIZE, cudaMemcpyDeviceToHost);
		
		//for logging - begin
		/*cudaMemcpy(hlog_indices, dlog_indices, 10*sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_positions, dlog_positions, 30*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_strengths, dlog_strengths, 30*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_c, dlog_c, 10*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_pos1, dlog_pos1, 30*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_pos2, dlog_pos2, 30*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_magnitude, dlog_magnitude, 10*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_springk, dlog_springk, 10*sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(hlog_length, dlog_length, 10*sizeof(float), cudaMemcpyDeviceToHost);
		PrintResults_StrengthKernel(fcuda, hlog_indices, hlog_positions, hlog_strengths, hlog_c,  hlog_pos1, hlog_pos2, hlog_magnitude,  hlog_springk, hlog_length);*/
		//for logging -end
		
		cudaFree(d_sources);
		cudaFree(d_strengths);
			
		//fclose(fcuda);
		
	}

}

//********************************************************** MEX  FUNCTION *****************************************************************
void mexFunction(int nlhs, mxArray *plhs[],   int nrhs, const mxArray *prhs[])
{    
	int devID;
	cudaDeviceProp props;
	cudaGetDevice(&devID);
	cudaGetDeviceProperties(&props, devID);

	//fmex = fopen("logmex.txt","a");
	
	double *p_sources 				= mxGetPr(prhs[0]);   
	double *p_num_stokeslets 		= mxGetPr(prhs[1]);   
	double *p_resting_length		= mxGetPr(prhs[2]);   
	double *p_vertical_length		= mxGetPr(prhs[3]); 
	double *p_diagonal_length		= mxGetPr(prhs[4]);   
	double *p_spring_constant	 	= mxGetPr(prhs[5]);
	double *p_num_towers	 		= mxGetPr(prhs[6]);
	
	num_sources 		= 		mxGetM(prhs[0]); 	// sources   
	          

	//host and pointer declarations
	float resting_length, vertical_length, diagonal_length,  spring_constant;
	float fnum_towers;
	float *h_sources, *h_strengths;
	
	//getting constant data from matlab pointers
	resting_length 		= *p_resting_length;		
	vertical_length 		= *p_vertical_length;
	diagonal_length		= *p_diagonal_length;
	spring_constant 	= *p_spring_constant;
	fnum_towers			= *p_num_towers;
	num_towers			= (int) fnum_towers;
	


	//allocate memory to host pointers
	ROWS 		= (int)(num_sources);
	COLS 		= 3;
	h_sources	= new float [ROWS * COLS  ];
	h_strengths	= new float [ROWS * COLS  ];
	memset(h_sources , 0, (ROWS*COLS)*sizeof(float));
	memset(h_strengths , 0, (ROWS*COLS)*sizeof(float));

		
	//copy data from matlab
	int t=0;		
	for (int ii=0;ii<COLS;ii++){				//no of columns
		for (int j=0;j<ROWS;j++){				//no of rows
			int i = col_array[ii];
			h_sources[ j*COLS + i]		= float(*(p_sources+t));
			t++;
		}
	}

	//print data recieved from  matlab
	//PrintFromMatlab(fmex, h_sources, resting_length, vertical_length, diagonal_length, spring_constant,  num_sources);

	ForceEval(h_sources, h_strengths, resting_length, vertical_length, diagonal_length, spring_constant);

	//print in host, what is being sent back to matlab
	//PrintToMatlab(fmex, h_strengths);	
	
	//paste data in matlab pointers	
	plhs[0] = mxCreateDoubleMatrix(num_sources, 3, mxREAL);
    double *ptrMatlab = mxGetPr(plhs[0]);		
	t=0;	
	for (int ii=0;ii<COLS;ii++){				//no of columns
		for (int j=0;j<ROWS;j++){				//no of rows
			int i = col_array[ii];
			*(ptrMatlab +t) = h_strengths[ j*COLS + i];
			t++;				
		}
	}

	//fprintf(fmex,"\n Returning from mex file succesfully ");	
	
	delete [] h_sources; 		// uninitialize
	delete [] h_strengths; 

	//fclose(fmex);
	
}


//********************************************************************Log Printing  Functions***********************************************************************************************************


void PrintResults_StrengthKernel(FILE *fp, int * hlog_indices, float * hlog_positions, float * hlog_strengths, float * hlog_c, float *  hlog_pos1, float *  hlog_pos2, float *   hlog_magnitude, float *  hlog_springk, float *  hlog_length){
	fprintf(fp, "\n  Strengths kernel is returning following values:  \n");
	fprintf(fp,"\nIndex,    					Position      		 c        			pos1  	  pos2			magnitude			springk				length				Strength : \n\n ");
	for(int i=0;i<10;i++)
		fprintf(fp," %d       %f %f %f       %f       %f %f %f      %f %f %f     %f       %f     %f           %f  %f  %f\n ", hlog_indices[i], hlog_positions[3*i], hlog_positions[3*i+1], hlog_positions[3*i+2], hlog_c[i], hlog_pos1[3*i], hlog_pos1[3*i+1], hlog_pos1[3*i+2],  hlog_pos2[3*i], hlog_pos2[3*i+1], hlog_pos2[3*i+2], hlog_magnitude[i], hlog_springk[i], hlog_length[i], hlog_strengths[3*i], hlog_strengths[3*i+1], hlog_strengths[3*i+2] );
		
}

void PrintFromMatlab(FILE *fp, float *h_sources, float resting_length, float vertical_length, float diagonal_length, float spring_constant,	int  num_sources){

	fprintf(fp,"\n finished copying and printing of what has been copied..Calling Kernels after this.\n");
	fprintf(fp,"\n variables:\n");
	fprintf(fp,"\n num_sources: = %d", num_sources);
	fprintf(fp,"\n resting_length: = %f", resting_length);
	fprintf(fp,"\n diagonal_length: = %f", diagonal_length);
	fprintf(fp,"\n vertical_length: = %f", vertical_length);
	fprintf(fp,"\n spring_constant: = %f", spring_constant);
	
	Print2dArray(fp, "\n strengths = \n", h_sources, ROWS,  COLS);	
	
	
}

void Print2dArray(FILE *fp, char *str,  float * ptr, int rows, int cols){
	fprintf(fp, "%s ", str);
	for(int i=0; i<rows*cols; i+=3)
		fprintf(fp, "%f  %f  %f\n", ptr[i], ptr[i+1], ptr[i+2]);

}

void PrintToMatlab(FILE *fp, float *h_strengths){
	fprintf(fp,"\n\n\nfinished kernel call and printing of what has been evaluated and being sent to matlab...\n\n");	
	Print2dArray(fp, "\n strengths = \n", h_strengths, ROWS,  COLS);	
}


