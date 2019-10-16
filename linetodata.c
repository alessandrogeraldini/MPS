// This module takes a line from a data file and converts it to an array of numbers (in practice, to a double pointer)
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Headerfiles/linetodata.h"

double *linetodata(char line[], int lenline, int *size)
{
	char testspace[2], numb[20];
	int i, lennumb, start_index, numb_index, space_index[lenline], words=0, index;
	double *numb_array, numb_double;
	/* Iterate over number of characters in the whole line
	 */
	for (i=0; i<lenline; i++)
	{
		testspace[0] = line[i];
		/* Compare each character to a white space. When a white space is found (or a 
		 * newline character), we break the string and extract the number from the 
		 * previous whitepace to the new one.
		 */
		if ( (strncmp(testspace, "\n", 1) == 0) || (strncmp(testspace, " ", 1) == 0) ) 
		{
			space_index[words] = i;
			words += 1;
		}
	}
	*size = words;
	numb_array = malloc(words*sizeof(double));
	for (index = 0;index < words; index++)
	{
		if (index == 0)
		{
			lennumb = space_index[index];
			start_index = 0;
		}
		else
		{
			lennumb = space_index[index] - space_index[index-1];
			start_index = space_index[index-1];
		}

	for (numb_index=0; numb_index < lennumb; numb_index++)
	{
		numb[numb_index] = line[start_index+numb_index];
	}
	for (numb_index=lennumb; numb_index < 20; numb_index++)
	{
		numb[numb_index] = ' ';
	}
 	sscanf(numb, "%lf", &numb_double);
	numb_array[index] = numb_double;
	}
return numb_array;
}
