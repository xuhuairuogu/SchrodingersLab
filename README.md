# HighNLSE
HighNLSE is a set of high order numerical solvers for the nonlinear Schrodinger equation using a pseudospectral method (Split step Fourier method).

The first function is \texttt{nlse_2.m}, which is your average second order solver.

The fourth, sixth and eighth order solvers have 2 versions: \texttt{nlse_n.m} and \texttt{nlse_n_NS.m}. The first one is much quicker but relies on subtraction, which might cause some numerical errors. The second one is slower but has no subtractions. It is quite slow especially for the eighth order.

The solvers also test for the energy error and plot it.

Note: this is a very incomplete work. If you do not know what you are doing and are not sure what results to expect, do not use it.
