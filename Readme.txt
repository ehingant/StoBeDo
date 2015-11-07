# StoBeDo
Simulation of the Stochastic Becker-Döring equations

This code is design to simulate the Stochastic Becker-Doring model.
It is a continuous time Markov chain whose generator and mathematical properties are
discribed in http://arxiv.org/abs/1412.5025
The algortihm is based on the Next reaction method presented in the paper:
    "Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many Channels
    Michael A. Gibson and and Jehoshua Bruck
    The Journal of Physical Chemistry A 2000 104 (9), 1876-1889
    DOI: 10.1021/jp993732q"
 
The code needs to be compiled with a C compiler such as GCC (http://gcc.gnu.org/) and linked with the GNU Scientific Library (http://www.gnu.org/software/gsl/) to works 

To compil just use in a terminal the command "gcc *.c -lgsl -lgslcblas -lm" in the directory with soure files.

Some functionality requier Gnulpot (http://www.gnuplot.info/) and ffmpeg (http://www.ffmpeg.org/).
 
The program is under licensed GNU GPLv3 (see a copy at http://www.gnu.org/licenses/gpl-3.0.en.html or in the directory).
It has been writting by E. Hingant (erwan 'at' mat.ufcg.edu.br) and R. Yvinec (romain.yvinec 'at' tours.inra.fr) between 2013 and 2015.

E.H would thanks Fondecyt Grant 3130318 (Chile) and CAPES/IMPA (Brazil) for their support.

