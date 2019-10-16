MPS_it: MPS_iterate.o MPS_single.o bilin_int.o linetodata.o tophat.o iondens.o newguess.o 
	gcc -O2 -Wall -ffast-math -o MPS MPS_iterate.o bilin_int.o linetodata.o tophat.o iondens.o newguess.o 
	gcc -O2 -Wall -ffast-math -o MPS_one MPS_single.o bilin_int.o linetodata.o tophat.o iondens.o newguess.o 
 
MPS_iterate.o: MPS_iterate.c 
	gcc -Wall -I/opt/local/include/gsl -c -o MPS_iterate.o MPS_iterate.c

MPS_single.o: MPS_single.c 
	gcc -Wall -I/opt/local/include/gsl -c -o MPS_single.o MPS_single.c

bilin_int.o: bilin_int.c
	gcc -Wall -I/opt/local/include/gsl -c -o bilin_int.o bilin_int.c

iondens.o: iondens.c
	gcc -Wall -I/opt/local/include/gsl -c -o iondens.o iondens.c

newguess.o: newguess.c
	gcc -Wall -I/opt/local/include/gsl -c -o newguess.o newguess.c

tophat.o: tophat.c
	gcc -Wall -I/opt/local/include/gsl -c -o tophat.o tophat.c

linetodata.o: linetodata.c
	gcc -Wall -I/opt/local/include/gsl -c -o linetodata.o linetodata.c
