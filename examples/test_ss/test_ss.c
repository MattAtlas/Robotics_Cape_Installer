/*******************************************************************************
* test_filters.c
*
* This demonstrates the use of the discrete time SISO filters. It sets up three
* filters, a complentary low & high pass filter along with an integrator.
* It varies a common input u from 0 to 1 through time and show the output of 
* each filter. It also displays the sum of the complementary high and low pass
* filters to demonstrate how they sum to 1
*******************************************************************************/

#include <robotics_cape.h>
#include <useful_includes.h>

#define SAMPLE_RATE 	50.0

int main(){

	printf("\n\nFirst lets start with a transfer function's polynomials.\n");
	float b_array[2] = {1, 5};
	vector_t b = create_vector_from_array(2, b_array);	
	float a_array[4] = {1, 5, 5, 5};
	vector_t a = create_vector_from_array(4, a_array);

	printf("\nb = \n");
	print_vector(b);
	printf("\na = \n");
	print_vector(a);

	printf("\nTo get this into state space form, use tf2ss.\n");
	CT_SS_filter_t CT_sys = tf2ss(b,a);

	printf("\nA = \n");
	print_matrix(CT_sys.A);
	printf("\nB = \n");
	print_matrix(CT_sys.B);
	printf("\nC = \n");
	print_matrix(CT_sys.C);


	printf("\nThis is of not much use with robots until put in discrete time.");
	printf("\nTo do this create a DT filter with the CT filter and dt as inputs.\n");
	DT_SS_filter_t DT_sys = create_DT_SS_filter(&CT_sys, 1/SAMPLE_RATE);

	printf("\nF = \n");
	print_matrix_sci_notation(DT_sys.F);
	printf("\nG = \n");
	print_matrix_sci_notation(DT_sys.G);
	printf("\nH = \n");
	print_matrix(DT_sys.H);


	printf("\nNow we want this filter to do something.  Let's make a control input.\n");
	vector_t u = create_vector(1);
	u.data[0] = 1;
	printf("\nu = \n");
	print_vector(u);

	printf("\nUsing the march_filter function, we give it the input and see how it reacts.\n");
	
	int counter;
	printf("  input u |");
	printf("  output  ");
	printf("\n");
	while(get_state() != EXITING){
	
		march_SS_filter(&DT_sys, u);
		
		printf("\r");
		printf("%7.2f   |", u.data[0]);
		printf("%7.2f   ", DT_sys.Xnew.data[0]*5);
		fflush(stdout);
		
		usleep(1000000/SAMPLE_RATE*2);
		
		// toggle u between 0 and 1 every 2 seconds
		counter++;
		if(counter >= SAMPLE_RATE){
			counter = 0;
			if(u.data[0]>0) u.data[0] = 0.0;
			else u.data[0] = 1.0;
		}
		
	}
	
	
	return 0;
}
