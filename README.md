# Magnetic Presheath Solver 

Download all files into a folder, then type `make` on the command line.

The input file _inputfile.txt_ has two values on two lines: alpha (in radians); = tau = T_i / T_e.

Set these to desired values: alpha = 0.05 and tau = 1 are good to start and are the default values after download.

Run the python scripts to produce phidata.txt and distfuncin.txt by typing `python phigenerator.txt` and `python Fgenerator.txt`.

Run the program by writing `./MPS`, which calculates the ion density profile and corrects the electrostatic potential contained in _phidata.txt_ until the quasineutrality equation _ion density = electron density_ is satisfied.

Writing `./MPS_one` just evaluates the ion density profile and the ion distributiion function at x=0 using the current _phidata.txt_.
