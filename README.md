# tsfcoll
tsfcoll

The authors are: Angelamaria Cardone, Dajana Conte, Beatrice Paternoster (Department of Mathematics, University of Salerno, Italy).

Some information about the authors:

Angelamaria Cardone webpage https://docenti.unisa.it/005020/home; email encardone@unisa.it
Dajana Conte webpage https://docenti.unisa.it/020280/home; email dajconte@unisa.it;
Beatrice Paternoster webpage https://docenti.unisa.it/000793/home; email beapat@unisa.it.

# Description


TSFCOLL solves an initial value problem for a fractional differential equation (FDE) by means of two step spline collocation methods.

 Two step spline collocation methods based are a generalization to FDEs of two step collocation methods for ODEs. This code implements methods introduced in [3], and further analyzed in [2]. A multivalue technique and a suitable graded mesh are used to obtain high order of convergence.  The starting procedure consists of the one step spline collocation methods  (compare [4] and https://github.com/cardange/tsfcoll) Please, cite this code as [1,2,3] if you need.

   [t,y,fval]=tsfcoll(f,b,gam,alpha,eta,r,N)
   integrates the initial value problem for the scalar FDE of order alpha>0
      D^alpha y(t) = f(t,y(t)), t in [0,b]
      y^(k)(0) = gam(k+1),      k=0,...,K-1
where K is the smallest integer grater than alpha, and D^alpha is the fractional derivative according to the Caputo's definition. f is a function handle corresponding to the vector field of the FDE for  the scalar time variable t and the state variable y.
The set of initial conditions gam is a vector of lenght K. 
eta(m) is equal to the vector of the collocation parameters. Constraint: 0 ≤ eta(1) ≤ · · · ≤eta(m).
r>=1 is grading exponent. For r=1 the mesh is uniform. For the optimal setting of r, compare [1,3].
N is the number of mesh points.

  IN OUTPUT
   t: vector of the graded mesh: t(i+1)=b(i/N)^r, i=0,...,N.
   y: vector of the approximate solution at the time steps t.
   fval is the number of the function evaluation. 

  REFERENCES

  [1] A. Cardone, D. Conte, B. Paternoster, (2023). A MATLAB Code for Fractional Differential Equations Based on Two-Step Spline Collocation Methods. In: Cardone, A., Donatelli, M., Durastante, F., 
Garrappa, R., Mazza, M., Popolizio, M. (eds) Fractional Differential Equations. INDAM 2021. Springer INdAM Series, vol. 50, 121–146. Springer, Singapore.  https://doi.org/10.1007/978-981-19-7716-9_8

[2] A. Cardone, D. Conte, B. Paternoster, Stability of two-step spline collocation methods for initial value problems for fractional differential equations,  Commun. Nonlinear Sci. Numer. Simul. 115, (2022), 106726.

[3] A. Cardone, D. Conte, B. Paternoster, Two-step collocation methods for fractional differential equations, Discrete Contin. Dyn. Syst. Ser. B 23 (2018), no. 7, 2709–2725. 

[4] Cardone A., Conte D., Paternoster B. (2021) A MATLAB Implementation of Spline Collocation Methods for Fractional Differential Equations. Lect. Notes Comput. Sci., vol 12949, 387-401.
