# HighNLSE
##IMPORTANT NOTE:
I tested this solved in MATLAB 2013a and 2015a and it had some issues with finite background waves in 2013a. This does not show in 2015a whatsoever, no idea what caused it.
This repo hasn't been maintained in a while, but I am working on it again and updating the solver.

##Introduction
HighNLSE is a set of high order numerical solvers for the nonlinear Schrodinger equation using the split-step Fourier method. A very simple explanation is available here: https://en.wikipedia.org/wiki/Split-step_method

##How to use it
To use this solver, you have to run the main script and modify the simulation parameters.

##Types of solvers:
###There are two classes of solvers:
1. Multi-product integrators: these are significantly quicker. See this paper for more info: http://arxiv.org/abs/0809.0914v2
2. Symplectic integrators: much slower than multi-product integrators but they do not contain subtraction, which might introduce an extra error. I highly recommend Multi-product integrators. See this paper for more on this method and the coefficients I used: http://journals.aps.org/pre/abstract/10.1103/PhysRevE.62.8746

###The available solvers are:

1. **T2**: Second-order splitting
2. **T4**: Fourth order splitting (multi-product integrators).
3. **T4_NS**: Fourth order splitting (symplectic integrators).
4. **T6**: Sixth order splitting (multi-product integrators).
5. **T6_NS**: Sixth order splitting (symplectic integrators).
6. **T8**: Eighth order splitting (multi-product integrators).
7. **T8_NS**: Eighth order splitting (symplectic integrators).

**Note**: the eighth order algorithms are actually too precise for double precision. They are useless on MATLAB, you might want to implement them in C, C++ or Fortran with quadruple precision

##Extra files
In addition, some extra files come with the solver and might not be of use to you. These were coded mainly for the current paper I'm working on and they might be of use to some people interested in verifying those results, or to other interested individuals. They are as follows:
1. ab: this compares the result to an Akhmediev breather.
2. energy: this computes the energy at every time step and plots the energy error dE. Useful for checking the order of the algorithm. You can verify this by plotting the integrated energy error fits to dE = a*dt^n where n is the order. This won't be a perfect fit, of course, and it won't really work for 8th order since you won't see its true precision in double.
3. a_plot: hard to explain, more details in file
4. recon: hard to explain, more details in file

**Note**: the last two files are very specific to my paper. Second one might be quite useful to you.

**Final note:** I coded this program in my sophomore undergraduate year as an electrical engineering major. I have done every effort I can to make sure that it works, but after all, my knowledge is very limited and it might have some mistakes. I'm just a student and I'm trying to learn. If you find anything, please let me know so I can fix it.

**This library/program is licenced under GNU GPL. The licence is available in the repo.**

Copyright 2015 Omar Ashour.
