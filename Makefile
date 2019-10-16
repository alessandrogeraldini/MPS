MPS_it: MPS_iterate.o MPS_single.o bilin_int.o linetodata.o tophat.o iondens.o newguess.o 
	gcc -O2 -Wall -ffast-math -o MPS ofiles/MPS_iterate.o ofiles/bilin_int.o ofiles/linetodata.o ofiles/tophat.o ofiles/iondens.o ofiles/newguess.o 
	gcc -O2 -Wall -ffast-math -o MPS_one ofiles/MPS_single.o ofiles/bilin_int.o ofiles/linetodata.o ofiles/tophat.o ofiles/iondens.o ofiles/newguess.o 
 
MPS_iterate.o: MPS_iterate.c 
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/MPS_iterate.o MPS_iterate.c

MPS_single.o: MPS_single.c 
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/MPS_single.o MPS_single.c

bilin_int.o: bilin_int.c
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/bilin_int.o bilin_int.c

iondens.o: iondens.c
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/iondens.o iondens.c

newguess.o: newguess.c
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/newguess.o newguess.c

tophat.o: tophat.c
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/tophat.o tophat.c

linetodata.o: linetodata.c
	gcc -Wall -I/opt/local/include/gsl -c -o ofiles/linetodata.o linetodata.c
