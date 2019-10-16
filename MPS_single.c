// Main function that calls the function to evaluate the ion density once, also evaluating the distributioin function at the Debye sheath entrance (x=0)
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "Headerfiles/iondens.h"

int main() 
{
clock_t begin_it = clock(); // Finds the start time of the computation
double tot_time;
iondens(1); // call ion density evaluation function
clock_t end_it = clock(); // finds the end time 
tot_time = (double) (end_it - begin_it) / CLOCKS_PER_SEC;
exit(0); 
}
