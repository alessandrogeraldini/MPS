// MODIFIED OCT 2019
// This is the main function of the iteration scheme to find the electrostatic potential in the magnetic presheath.
// It calls the function to evaluate the ion density, then calls the function to evaluate the new guess.
// If the ion density is sufficiently close to the electron density, the function to evaluate the ion density is called one last time and the distribution function at x=0 (the Debye sheath entrance) is computed. 
// All functions called by this main function are in separate files.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Headerfiles/iondens.h"
#include "Headerfiles/newguess.h"

int main() 
{
clock_t begin_it = clock(); // Finds the start time of the computation
int convergence = 0, problem = 0; 
/* convergence = 0 tells newguess.c to set smoothing for the next iteration, convergence = 1 tells it to not smooth, and convergence = 2 means that the convergence criterion has been satisfied and the iteration is finished */
int N_iterations = 1000, N;
/* Number of maximum iterations is set to some large number e.g. 100 but iterations should converge in about 20. If they don't, then there is a problem. */
double tot_time;
srand ( time ( NULL));
for (N=0;N<N_iterations;N++)
{	
	printf("Iteration number is %d\n", N);
	problem = iondens(0); //The argument of iondens set to zero makes the module not compute the ion distribution function
	convergence = newguess(convergence, problem);
	if (convergence == 2)
	{	
		iondens(1); //The argument of iondens set to one makes the module compute the ion distribution function at x=0
		clock_t end_it = clock(); // finds end time of last iteration
		tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
		printf("At %dth iteration converged successfully in %f seconds\n", N, tot_time);
		exit(0); 
	} 
}
if (N_iterations == 0)
{	
	problem = iondens(1);
}
printf("No convergence after %d iterations\n", N_iterations);
exit(0);
}
