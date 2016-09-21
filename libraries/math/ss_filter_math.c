/*******************************************************************************
* filter_math.c
*
* Matt Atlas 2016
*******************************************************************************/

#include "../robotics_cape.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DEBUG
#define PI 3.14159265

/*******************************************************************************
* CT_SS_filter_t createCTSSfilter(int states, int inputs, int outputs)
*
* Sets up all matrices for system
*******************************************************************************/
CT_SS_filter_t create_CT_SS_filter(int states, int inputs, int outputs){

	CT_SS_filter_t filter;

	if((states < 1) | (inputs < 1) | (outputs < 1)){
		printf("Error: Invalid values.\n");
		return filter;
	}
	
	filter.states  = states;
	filter.inputs  = inputs;
	filter.outputs = outputs;
	
	filter.A = create_square_matrix(filter.states);
	filter.B = create_matrix(filter.states,filter.inputs);
	filter.C = create_matrix(filter.outputs,filter.states);
	filter.D = create_matrix(filter.outputs,filter.inputs);
	
	filter.Xold = create_vector(filter.states);
	filter.Xnew = create_vector(filter.states);
	filter.Y    = create_vector(filter.outputs);

	
	filter.saturation_en   = 0;
	filter.saturation_flag = 0;
	filter.saturation_high = create_vector(filter.inputs);
	filter.saturation_low  = create_vector(filter.inputs);
	
	filter.K = create_matrix(filter.inputs,filter.states);
	filter.L = create_matrix(filter.states,filter.outputs);
	
	return filter;
}


/*******************************************************************************
* DT_SS_filter_t createDTSSfilter(CT_SS_filter_t* CT_sys, float dt)
*
* Sets up all matrices for system
*******************************************************************************/
DT_SS_filter_t create_DT_SS_filter(CT_SS_filter_t* CT_sys, float dt){
	DT_SS_filter_t filter;
	// error handling
	filter.dt = dt;

	filter.states  = CT_sys->A.rows;
	filter.inputs  = CT_sys->B.cols;
	filter.outputs = CT_sys->C.rows;
	
	filter.Xold = create_vector(filter.states);
	filter.Xnew = create_vector(filter.states);
	filter.Y    = create_vector(filter.outputs);
	
	filter.F = create_square_matrix(filter.states);
	filter.G = create_matrix(filter.states,filter.inputs);
	filter.H = create_matrix(filter.outputs,filter.states);
	filter.D = create_matrix(filter.outputs,filter.inputs);

	// DT model of system	
	filter.F = C2D_A2F(CT_sys->A, filter.dt);
	filter.G = C2D_B2G(CT_sys->A, CT_sys->B, filter.dt);
	filter.H = CT_sys->C;
	
	filter.saturation_en   = 0;
	filter.saturation_flag = 0;
	filter.saturation_high = create_vector(filter.inputs);
	filter.saturation_low  = create_vector(filter.inputs);
	
	filter.K = create_matrix(filter.inputs,filter.states);
	filter.L = create_matrix(filter.states,filter.outputs);
	
	return filter;
}



/*******************************************************************************
* int tf2ss(vector_t b, vector_t a, CT_SS_filter_t* CT_sys)
* Takes transfer function polynomials and converts into state space
* continuous time controller canonical form.
*******************************************************************************/
CT_SS_filter_t tf2ss(vector_t b, vector_t a){
	int i;

	CT_SS_filter_t CT_sys = create_CT_SS_filter(a.len-1, 1, 1);
	// Make sure it's proper or semiproper.
	if (b.len > a.len){
		printf("Error: Improper transfer function\n");
		return CT_sys;
	}
	if(a.data[0] == 0.0){
		printf("Error: Leading coefficient on 'a' can't be zero\n");
		return CT_sys;
	}
	// Make a new b_temp vector the same length as a,
	// and allocate memory with zero coefficients.
	vector_t b_temp = create_vector(a.len);

	// Put the values of b in b_temp starting at the end and moving back.
	// This is to make leading coeffs zero if a > b.
	for(i=1;i<=b.len;i++){
		b_temp.data[a.len - i] = b.data[b.len - i];
	}
	// If leading coefficient on a is not 1, normalize transfer function.
	if(a.data[0] != 1){
		float coeff = a.data[0];
		for(i=0;i<=a.len;i++){	
			a.data[i] = a.data[i] / coeff;	
			b_temp.data[i] = b_temp.data[i] / coeff;
		}
	}
	// For a semiproper transfer function, we need to subtract off terms from C
	// and put in the D term to CT_sys.
	if(b_temp.data[0] != 0){
		for(i=0;i<a.len;i++){
			b_temp.data[i] = b_temp.data[i] - a.data[i+1]*b_temp.data[0];
		}
		CT_sys.D.data[0][0] = b_temp.data[0];
	}

	// Fill in A matrix
	for(i=0;i<a.len-2;i++){
		CT_sys.A.data[i+1][i] = 1;				// fill in lower identity
	}
	for(i=0;i<a.len-1;i++){
		CT_sys.A.data[0][i] = -a.data[i+1];	// fill in top row
	}
	// Fill in B matrix
	CT_sys.B.data[0][0] = 1;
	// Fill in C matrix
	for(i=0;i<a.len-1;i++){
		CT_sys.C.data[0][i] = b_temp.data[i+1];
	}
	return CT_sys;
}

/*******************************************************************************
* 
*
* 
*******************************************************************************/
// int saturate(DT_SS_filter_t* sys, matrix_t* input){
	
	// for (int i=0;i<sys->inputs;i++){
		// if (&input[i] > sys->saturation_high[i]){
			// &input[i] = sys->saturation_high[i];
			
			// sys->saturation_flag = 1;
		// }
		// else if(&input[i] < sys->saturation_low[i]){
				// &input[i] = sys->saturation_low[i];
			
			// sys->saturation_flag = 1;
		// }
		// else{
			// sys->saturation_flag = 0;
		// }
	// }
	
	// return 0;
// }

/*******************************************************************************
* int march_filter(DT_SS_filter_t* DT_sys, vector_t input)
*
* 
*******************************************************************************/
int march_SS_filter(DT_SS_filter_t* DT_sys, vector_t input){
	
	vector_t FX;	// temporary vectors
	vector_t Gu;
	int i;
	
	DT_sys->Xold = DT_sys->Xnew;
	
	if(input.len != DT_sys->G.cols){
		printf("Error: input vector size mismatch");
		return -1;
	}
	
	if(DT_sys->saturation_en == 1){
		for(i=0;i<input.len;i++){
			if(input.data[i] < DT_sys->saturation_low.data[i]){
				input.data[i] = DT_sys->saturation_low.data[i];
				DT_sys->saturation_flag = 1;
			}
			if(input.data[i] > DT_sys->saturation_high.data[i]){
				input.data[i] = DT_sys->saturation_high.data[i];
				DT_sys->saturation_flag = 1;
			}
		}
	}
	
	FX = matrix_times_col_vec(DT_sys->F, DT_sys->Xold);
	Gu = matrix_times_col_vec(DT_sys->G, input);
	
	for(i=0;i<input.len;i++){
		DT_sys->Xnew.data[i] = FX.data[i] + Gu.data[i];
	}

	DT_sys->Y = matrix_times_col_vec(DT_sys->H, DT_sys->Xnew);
	#ifdef DEBUG
	// output y is estimate of sensor data
	//printf("y = %f\n",DT_sys->Xold.data[2][0]*DT_sys->H.data[0][2]);
	#endif

	destroy_vector(&FX);
	destroy_vector(&Gu);
	
	return 0;
}



/*******************************************************************************
* matrix_t C2D_A2F(matrix_t A, float h)
*
* Using CT A matrix and time step, get DT F matrix
*******************************************************************************/
matrix_t C2D_A2F(matrix_t A, float h){
	int m = A.rows;
	int i,j,k,N;
	matrix_t sumold = create_square_matrix(m);
	matrix_t sumnew = create_square_matrix(m);
	matrix_t result = create_square_matrix(m);
	float sum;

	// initialize identity matrix for sumold and result
	for(i=0;i<m;i++){		
		sumold.data[i][i] = 1.0;	// A^0 first element of sum
		result.data[i][i] = 1.0;
	}
	// N = order of exponential expansion sum
	for(N=1;N<5;N++){

		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				sum = 0; 	// initialize sum for next loop
				for(k=0;k<m;k++){
					// do the matrix multiplication
					sum += sumold.data[i][k]*A.data[k][j];
				}
				// save mult sum to new location
				sumnew.data[i][j] = h*sum/N;
			}
		}	

		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				//matrix exponential expansion sum over N
				result.data[i][j] += sumnew.data[i][j];				
				sumold.data[i][j] = sumnew.data[i][j];	
			}
		}
	}
	return result;
}

	
/*******************************************************************************
* matrix_t C2D_B2G(matrix_t A, matrix_t B, float h)
*
* Using CT A and B matrices and time step, get DT G matrix
*******************************************************************************/
matrix_t C2D_B2G(matrix_t A, matrix_t B, float h){
	int i,j,k,N;
	int m = A.rows;
	matrix_t sumold = create_square_matrix(m);
	matrix_t sumnew = create_square_matrix(m);
	matrix_t result = create_square_matrix(m);
	float sum;
	matrix_t G = create_matrix(B.rows,B.cols);

	// initialize identity matrix
	for(i=0;i<m;i++){		
		sumold.data[i][i] = 1;	// A^0 first element of sum
		result.data[i][i] = 1;
	}
	// N = order of exponential expansion sum
	for(N=1;N<5;N++){				
		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				sum = 0; 			// initialize sum for next loop
				for(k=0;k<m;k++){
					// do the matrix multiplication
					sum += sumold.data[i][k]*A.data[k][j];
				}
				// save mult sum to new location
				sumnew.data[i][j] = h*sum/(N+1);// factorial starts 1 higher than in F
			}
		}	
		for(i=0;i<m;i++){
			for(j=0;j<m;j++){
				//matrix exponential expansion sum over N
				result.data[i][j] += sumnew.data[i][j];				
				sumold.data[i][j] = sumnew.data[i][j];	
			}
		}
	}
	for(i=0;i<m;i++){
		sum = 0;
		for(j=0;j<B.cols;j++){
			for(k=0;k<B.rows;k++){
			sum += result.data[i][k]*B.data[k][j];
			}
		G.data[i][j] = h*sum;
		}
	}
	return G;
}

