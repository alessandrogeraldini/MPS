// LAST MODIFIED OCT 2019
// This code calculates the next electrostatic potential guess in the iteration to obtain the self-consistent magnetic presheath electrostatic potential profile 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "Headerfiles/linetodata.h"
#include "Headerfiles/newguess.h"
#include "Headerfiles/bilin_int.h"

int newguess(int convergence, int problem)
{//DECLARATIONS
	clock_t begin = clock(); // Finds the start time of the computation
	int debug = 0;
	int row, col, i, ii, sizemumu, sizeUU, j, L1=0, n;
	char line[200000];
	double *storevals, *gg, *xx, *ni, **FF, *mumu, *UU, k32, k32denom, k32num, k1DK, *Tevect, *k1GKvect; 
	double *phipg, *phippg, *phip, *phipp;;
	double k1DKnum1old, densinf, fluxinf, k1GK, **FFprime, k32denom1, k1DKnum1, k1DKnum, *newphi, k32denom1old, fluxinf1old, fluxinf1;
	double densinf1old, densinf1;
	double CC, *flux, *ne, k1DKthres = 0.2, Fprime, Fprimeold, F, Fold, Te, Teor1, fluxinfintgrdold, fluxinfintgrd;
	double phi1, alpha, weight, U, u, du, *phi;
	double smoothweight = 2.0/3.0;
	double sig=0.0, siglimit, siglimit1, dev, devbig;
	double random_value;


	FILE *distfile, *densfile, *Umufile, *fout2, *input, *fp, *fpotold;
if (problem == 1) 
{// If a problem was found evaluating the ion density, set convergence back to zero for smoothing to re-start
	convergence = 0; 
}
if ((fp = fopen("phidata.txt", "r")) == NULL)
{// Check for presence of file
	printf("Cannot open %s\n", "phidata.txt");
	exit(EXIT_FAILURE); }
// The while loop counts the lines in the file to determine the size of the arrays to be created 
row=0;
while(fgets(line, 20000, fp) != NULL)
{ // Count the number of rows in the distfuncin.txt file
	row += 1; 
}
n=row;
row=0;
// Below we initialize all arrays that contain functions of position x with the correct size n 
ne = calloc(n,sizeof(double)); // phi now has correct size
ni = calloc(n,sizeof(double)); // phi now has correct size
phi = calloc(n,sizeof(double)); // phi now has correct size
phip = calloc(n,sizeof(double)); // phi now has correct size
phipg = calloc(n,sizeof(double)); // phi now has correct size
phipp = calloc(n,sizeof(double)); // phi now has correct size
phippg = calloc(n,sizeof(double)); // phi now has correct size
gg = calloc(n,sizeof(double)); // gg = sqrt(x) has correct size
xx = calloc(n,sizeof(double)); // xx = x has correct size
newphi = calloc(n,sizeof(double)); // same as above 
flux = calloc(n,sizeof(double)); // same as above 
row = 0;
rewind(fp);
while (fgets(line, 2000, fp) != NULL)
{	//Fill phi and gg arrays
	storevals = linetodata(line, strlen(line), &col);
	gg[row] = *storevals;
	phi[row] = *(storevals+1);
	xx[row] = gg[row]*gg[row];
	if (debug == 1)
	{	
		printf("at row = %d, phi = %f, gg = %f, xx = %f\n", row, phi[row], gg[row], xx[row]); 
	}
	row += 1; 
}
fclose(fp);
if ((input = fopen("inputfile.txt", "r")) == NULL)
{ // Check for presence of file
       	printf("Cannot open %s\n", "inputfile.txt");
        exit(1); }
i=0;
while (fgets(line, 20, input) != NULL)
{	
        storevals = linetodata(line, strlen(line), &col);
	if (i==0)
	{	alpha = *storevals; }
	else if (i==1)
	{	Te = *storevals; }
        i += 1; }
fclose(input);
if (Te>1.001)
{	
	Teor1 = Te;
}
else 
{	
	Teor1 = 1.0; 
}
siglimit1 = 0.03;
siglimit = 0.007;
if (Te>2.1 || alpha < 0.01)
{
	siglimit1 = 0.02;
	siglimit = 0.012;
}
// Assign value to weight depending on whether convergence flag 

random_value = (double) rand()/RAND_MAX;//float in range 0 to 1
printf("The random number between 0 and 1 is %f\n", random_value);
if (convergence == 0)
{	
	weight = 0.8*random_value/pow(Teor1, 1.0); 
}
else
{	
	weight = 0.4*random_value/pow(Teor1, 1.0); 
}
if ((Umufile = fopen("Umufile.txt", "r")) == NULL)
{
	printf("cannot open file %s\n", "Umufile.txt");
	exit(-1);
}
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and assigning the values of the first line to mumu, second to UU  */
row=0;
while (fgets(line, 200000, Umufile) != NULL)
{	storevals = linetodata(line, strlen(line), &col);
	if (row == 0)
	{	mumu = storevals;
		sizemumu = col; 
	}
	else
	{	UU = storevals;
		sizeUU = col; 
	}
	row += 1; 
}
row=0;
fclose(Umufile);
/* This part of the code extracts the distribution function from a file. We import the distribution function into a 2 dimensional array, F(mu,U). */
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1); 
}
while(fgets(line, 20000, distfile) != NULL)
{ // Count the number of rows in the distfuncin.txt file
	row += 1; }
// Allocate the right amount of memory to the distribution function pointer
FF = malloc(row*sizeof(double));
// The number of rows is also the the size of the array in UU
row = 0; // Set number of rows counter to zero again
fclose(distfile); // Close file
if ((distfile = fopen("distfuncin.txt", "r")) == NULL)
{	// Check file exists and can be opened
	printf("cannot open file %s\n", "distfuncin.txt");
	exit(-1); }
while (fgets(line, 200000, distfile) != NULL)
{ /* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to FF and assign each number to a different column of FF (for a fixed row). The  while loop then jumps to a new line and repeats the process. */
	storevals = linetodata(line, strlen(line), &col);
	FF[row] = storevals;
	row +=1; 
}
fclose(distfile);
row = 0;
if ((densfile = fopen("PostProcessing/niout.txt", "r")) == NULL)
{	
	printf("Cannot open file %s\n", "niout.txt"); 
	exit(-1); 
}
while (fgets(line, 20000, densfile) != NULL)
{
/* While loop dedicated to reading each line of the file, extracting the data (ie the numbers in the line) and counting the columns. Once the columns are counted we allocate memory to FF and assign each number to a different column of FF (for a fixed row). The  while loop then jumps to a new line and repeats the process. */
	storevals = linetodata(line, strlen(line), &col);
	ni[row] = *(storevals+1);
	//printf("ni[%d] = %f\n", row, ni[row]);
	row += 1; 
}
fclose(densfile);
L1 = row;
printf("L1 = %d\n", L1);
if ((fpotold = fopen("PostProcessing/phidataold.txt", "w")) == NULL)
{
	printf("Cannot open phidataold.txt");
	exit(EXIT_FAILURE);
}
for (i=0; i<n; i++)
{
	fprintf(fpotold, "%f %f\n", gg[i], phi[i]);
}
fclose(fpotold); 
// INTEGRALS OF DISTRIBUTION FUNCTION AT INFINITY
/* In the section below I will evaluate the integrals that are necessary at the Chodura sheath entrance. They will be used in the potential iteration to extrapolate the electrostatic potential at large values of x which inevitably the code cannot solve for (boundary condition at infinity). */
FILE *fpk1;
if ((fpk1 = fopen("k1gyrokinetic.txt", "r")) == NULL)
{	// Check for presence of file
	printf("Cannot open %s\n", "k1gyrokinetic.txt");
	exit(EXIT_FAILURE); 
}
/* The while loop counts the lines in the file to determine the size of the arrays to be created */
row=0;
while(fgets(line, 20000, fpk1) != NULL)
{ // Count the number of rows in the k1gyrokinetic.txt file
	row += 1; 
}
k1GKvect = malloc(row*sizeof(double));
Tevect = malloc(row*sizeof(double));
rewind(fpk1);
row=0;
while (fgets(line, 2000, fpk1) != NULL)
{
	storevals = linetodata(line, strlen(line), &col);
	Tevect[row] = *storevals;
	k1GKvect[row] = *(storevals+1);
	//printf("i=%d, Tevect=%f and k1GKvect=%f\n", i, Tevect[i], k1GKvect[i]);
	row += 1;
}
fclose(fpk1); // close file
printf("Tevect = %f\n", Tevect[0]);
i=0;
while (Tevect[i] < Te)
{
//	printf("Tevect = %f\n", Tevect[i]);
	i+=1;
}
k1GK = ( k1GKvect[i]*(-Tevect[i-1] + Te) + k1GKvect[i-1]*(Tevect[i] - Te) )/(Tevect[i] - Tevect[i-1]);
//printf("k1GK = %f\n", k1GK);

// Note: first index of FF and FFprime is mu, second one is U

FFprime = malloc(sizemumu*sizeof(double));
//printf("FF[30][29] = %f and FF[29][30] =%f \n", FF[30][29], FF[29][30]);


/* The loop below evaluates the derivative of F with respect to total energy U at U = mu
 */
for (i=0; i<sizemumu; i++)
{
	FFprime[i] = malloc(sizeUU*sizeof(double));
	for (j=0; j < sizeUU; j++)
	{
		if ((j==sizeUU-1) || (i==sizemumu-1))
		{
			FFprime[i][j] = 0.0;	
		}
		else
		{
			FFprime[i][j] = (FF[i][j+1] - FF[i][j])/(UU[j+1] - UU[j]); 
		}
	//FFprime[i] = (1/pow(M_PI, 1.5))*exp(-mumu[i]);
		//printf("FFprime = %f\n", FFprime[i][j]);
	}
}

k32num = 0.0;
for (j=0; j<sizemumu; j++)
{
	if (j!=0)
	{
		k32num += 0.5*0.5*(16.0*M_PI/3)*(FFprime[j][0]+FFprime[j-1][0])*(mumu[j] - mumu[j-1]);
		//printf("k32num = %f\n", k32num);
	}
}
// Numerator ok

k1DKnum1old = 0.0;
k32denom1 = 0.0;
k1DKnum1 = 0.0;
k32denom1old = 0.0;
k32denom = k1DKnum = 0.0;
fluxinf1old = fluxinf1 = fluxinf = 0.0;
densinf1old = densinf1 = densinf = 0.0;
du = 0.01;
for (i=0; i<sizemumu; i++)
{
	if (i!=0)
	{
	k32denom1old = k32denom1;
	k32denom1 = 0.0;
	k1DKnum1old = k1DKnum1;
	k1DKnum1 = 0.0;
	fluxinf1old = fluxinf1;
	fluxinf1 = 0.0;
	densinf1old = densinf1;
	densinf1 = 0.0;
	}
	Fprimeold = Fprime = 0.0;
	Fold = F = 0.0;
	fluxinfintgrdold = fluxinfintgrd = 0.0;
	for (j=0; j< (int) 500*sqrt(Teor1); j++)
	{
		u = du*j;
		fluxinfintgrdold = fluxinfintgrd;
		Fprimeold = Fprime;
		Fold = F;
		U = pow(u, 2.0) + mumu[i];
		Fprime = bilin_interp(mumu[i], U-mumu[i], FFprime, mumu, UU, sizemumu, sizeUU, -1, -1);
		F = bilin_interp(mumu[i], U-mumu[i], FF, mumu, UU, sizemumu, sizeUU, -1, -1);
		fluxinfintgrd = F*u*alpha;
		k32denom1 += 0.5*(Fprimeold + Fprime)*du;
		k1DKnum1 += 0.5*(Fprimeold + Fprime)*du;
		fluxinf1 += 0.5*(fluxinfintgrd + fluxinfintgrdold)*du;
		densinf1 += 0.5*(F + Fold)*du;
		//printf("k32denom1 = %f\n", k32denom1);
	}
	k32denom1 *= mumu[i];
	if (i!=0)
	{
	k32denom += 0.5*(mumu[i]-mumu[i-1])*(k32denom1 + k32denom1old);
	k1DKnum += 0.5*(mumu[i]-mumu[i-1])*(k1DKnum1 + k1DKnum1old);
	fluxinf += 0.5*(mumu[i]-mumu[i-1])*(fluxinf1 + fluxinf1old);
	densinf += 0.5*(mumu[i]-mumu[i-1])*(densinf1 + densinf1old);
	}
}


k32denom = 1 + 4.0*M_PI*k32denom;
k1DKnum = 1 - 4.0*M_PI*k1DKnum*Te;

k32 = k32num/k32denom;
k1DK = 2.0*k1DKnum/(Te*k32denom);
k1DK *= Te;
densinf *= (4*M_PI);
fluxinf *= (4*M_PI); 

printf("k32num = %f\n", k32num);
printf("k32 = %f\n", k32);
printf("k1DK = %f, k1DK threshold is %f\n", k1DK, k1DKthres);
printf("Te = %f\n", Te);
printf("fluxinf = %f\n", fluxinf);

printf("Integrals of distribution function at infinity DONE\n");

sig = 0.0;
devbig = dev = 0.0;
// The part that calculates the new electrostatic potential
for (i=0; i<L1; i++)
{
	ne[i] = exp(phi[i]/Te);
	dev = fabs(ni[i]/ne[i] - 1.0);
	if (dev > devbig)
	{	devbig = dev; }
	sig += pow((ni[i]/ne[i] - 1.0), 2.0);
	if (i==L1-1)
	{
		//printf("ni= %f, weight = %f, phi[i] = %f\n", ni[i], weight, phi[i]);
		//newphi[i] = phi[i] - weight*(ne[i] - ni[i]);
		if (problem == 0)
		{
			newphi[i] = Te*log(ne[i] + weight*(ni[i] - ne[i]));
		}
		else
		{
			newphi[i] = 0.0;
		}
		//newphi[i] = Te*log(weight*ni[i] + (1-weight)*ne[i]);
		if (fabs(k1DK) < k1DKthres)
		{
			CC = sqrt(20.0/k32)/pow(-0.5*newphi[i], 0.25) - gg[i]*gg[i];
			printf("C_3/2 = %f, newphi[i] = %f\n", CC, newphi[i]);
		}
		else if (fabs(k1DK) > k1DKthres)
		{
			phi1 = newphi[i]*exp(-k1GK*xx[i]);
		}
	}
	else 
	{
		// This is where the iteration happens
		//newphi[i] = phi[i] - weight*(ne[i] - ni[i]);
		if (problem == 0)
		{
			newphi[i] = Te*log(ne[i] + weight*(ni[i] - ne[i]));
		}
		else
		{
			newphi[i] = 0.0; 
		}
		//newphi[i] = Te*log(weight*ni[i] + (1-weight)*ne[i]); (exactly the same as the above)
		//printf("ni= %f\n", ni[i]);
		//printf("newphi[i] = %f\n", newphi[i]);
	}
}
for (i=L1; i<n; i++)
{
	if (fabs(k1DK) < k1DKthres)
	{
		//printf("%f = newphi[i]\n", newphi[i]);
		newphi[i] = -2.0*(400.0/pow(k32, 2.0))/pow(gg[i]*gg[i] + CC, 4.0);
	}
	else if (fabs(k1DK) > k1DKthres)
	{
		newphi[i] = phi1*exp(k1GK*xx[i]); // k1GK should be negative
		//printf("phi = %f, newphi[i] = %f\n", phi[i], newphi[i]);
	}
	ne[i] = exp(phi[i]/Te);
	flux[i] = fluxinf/ni[i];
}
sig /= L1;
sig = sqrt(sig);

if (sig < siglimit && devbig < 5.0*siglimit) // || devbig < 0.04)
{
	convergence = 2;
	printf("convergence = %d, sig = %f\n", convergence, sig);
}
else
{
	// SMOOTHING OF PHI
	// Now we take derivatives of phi
	for (i=0; i<n; i++)
	{	
		if (i != n-1)
		{	phip[i] = (newphi[i+1] - newphi[i])/(gg[i+1]*gg[i+1] - gg[i]*gg[i]);
			phipg[i] = (newphi[i+1] - newphi[i])/(gg[i+1] - gg[i]);
			if (i!=0)
			{	phipp[i] = 2.0*(phip[i] - phip[i-1])/(gg[i+1]*gg[i+1] - gg[i-1]*gg[i-1]);
				phippg[i] = 2.0*(phipg[i] - phipg[i-1])/(gg[i+1] - gg[i-1]); }	
			if (i==1)
			{	phipp[0] = phipp[1];
				phippg[0] = phippg[1]; } }
		else
		{	
			phip[i] = phip[i-1];
			phipg[i] = phipg[i-1];
			phipp[i] = phipp[i-1];
			phippg[i] = phippg[i-1];
		} 
	}
	/* Smooth second derivative of phi wrt sqrt(x)
	*/
	if (convergence == 0)
	{ 	
		for (i=0; i<L1; i++)
		{	newphi[i] = phi[L1];
			phipg[i] = phipg[L1];
			if (i!=0)
			{	
				phippg[i] = ((1.0-smoothweight)*phippg[i-1]*(gg[i+1] - gg[i]) + smoothweight*phippg[i]*(gg[i+1] - gg[i-1]) + (1.0-smoothweight)*phippg[i+1]*(gg[i] - gg[i-1]))/((gg[i+1] - gg[i-1]));
				//phipp[i] = ((1.0-smoothweight)*phipp[i-1]*(xx[i+1] - xx[i]) + smoothweight*phipp[i]*(xx[i+1] - xx[i-1]) + (1.0-smoothweight)*phipp[i+1]*(xx[i] - xx[i-1]))/((xx[i+1] - xx[i-1])); 
			} 
		}
		// Re-integrate smoothed second derivative to first
		for (i=L1-1; i>=0; i--)
		{	for (ii=i; ii>=0; ii--)
			{	phipg[ii] -= (gg[i+1] - gg[i])*(phippg[i+1] + phippg[i])*0.5; 
			} 
		}  
				//phip[ii] -= (xx[i+1] - xx[i])*(phipp[i+1] + phipp[i])*0.5; 
		// Re-integrate first derivative to phi
		for (i=L1-1; i>=0; i--)
		{	for (ii=i; ii>=0; ii--)
			{	
				newphi[ii] -= (gg[i+1] - gg[i])*(phipg[i+1] + phipg[i])*0.5; 
			} 
		} 
	}
				//newphi[ii] -= (xx[i+1] - xx[i])*(phip[i+1] + phip[i])*0.5; 
	if ((fout2 = fopen("phidata.txt", "w")) == NULL)
	{
		printf("Cannot open phidata.txt");
		exit(EXIT_FAILURE);
	}
	for (i=0; i<n; i++)
	{	
		fprintf(fout2, "%f %f\n", gg[i], newphi[i]);
	}
	fclose(fout2);
	printf("convergence = %d, sig = %f, devbig = %f\n", convergence, sig, devbig);
	if (convergence == 0 && ( sig < siglimit1 ) ) //|| devbig < 0.2) )
	{
		convergence = 1; 
		printf("now don't smooth = %d, sig = %f\n", convergence, sig);
	}
	else if ((convergence == 1) && (sig > siglimit1))
	{
		convergence = 0;
	}
}
free(ne); free(ni); free(phi); free(phip); free(phipg); free(phipp); free(phippg); free(gg); free(xx); free(newphi); free(flux);
free(k1GKvect); free(Tevect);
for (row=0;row<sizemumu;row++)
{
	free(FF[row]);
	free(FFprime[row]);
}
free(FF); free(FFprime);
clock_t end = clock(); // finds the end time of the computation
double jobtime  = (double)(end - begin) / CLOCKS_PER_SEC;
printf("Module to obtain the new guess for the electrostatic potential ran in %f\n", jobtime);
return convergence;
}
