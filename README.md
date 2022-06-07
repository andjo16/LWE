By Andreas Holm Jørgensen  (s202743)
Used for Master Thesis project within the subject of Post Quantum Cryptography  
for completing a degree in Master of Science in Engineering (Computer Science) at Technical University of Denmark (DTU), June 2022

This folder contains the python scripts and result-files used for documenting the
empiricical decryption error rate for different set of parameters for the LWE encryption scheme.

Python files
 - `LWE64.py`: an implementation of the LWE encryption scheme
 - `FFPNG.py`: Used to find the impirical success probability of forcing decryption errors given we have an LWE oracle and therefore knows error-matrix E
 - `parameter instantiation.py`: used to find suitable sets of parameters as described in the report. It also uses LWE64.py and FFPNG.py 
 to estimate empirical rate of decryption errors
 
 Other files
 - `combinedResultsFinal.npy`: The file containing the data giving results presented in the report. 
 The file contains a numpy matrix of all the results, and us automatically loaded by `parameter instations.py`
 - `combinedLatexFinal.txt`: The latex-formatting of the data, in the form presented in the report.
 
The script runs under python 3. It can be run with the command  
`python "parameter instantiation.py"`

When the above is run, the file will recalculate the all the sets of parameters, load the data in `combinedResultsFinal.npy` 
and merge the data, before printing the results as a latex table.  
The script will generate two folders `results/` and `results2/`. These folders will only contain temporary backup files.

To generate empirical results anew, set the variables `findErrorProb`=`True` and `findForcedErrorProb`=`True`.  
However be aware that the script is expected to run for many hours if this is done.

 
 