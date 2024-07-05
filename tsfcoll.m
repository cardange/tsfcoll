function [t,y,fval]=tsfcoll(f,b,gam,alpha,eta,r,N)
%%  
%  TSFCOLL solves an initial value problem for a fractional differential equation 
%  (FDE) by means of two step spline collocation methods.
%
%  Two step spline collocation methods based are a generalization to FDEs 
%  of two step collocation methods for ODEs. This code implements methods 
%  introduced in [3], and further analyzed in [2]. A multivalue technique 
%  and a suitable graded mesh are used to obtain high order of convergence.
%  The starting procedure consists of the one step spline collocation methods
%  (compare [4] and https://github.com/cardange/tsfcoll)
%  Please, cite this code as [1,2,3] if you need.
%
%   [t,y,fval]=tsfcoll(f,b,gam,alpha,eta,r,N)
%   integrates the initial value problem for the scalar FDE of order alpha>0
%      D^alpha y(t) = f(t,y(t)), t in [0,b]
%      y^(k)(0) = gam(k+1),      k=0,...,K-1
%   where K is the smallest integer grater than alpha, and D^alpha is the
%   fractional derivative according to the Caputo's definition. f is a
%   function handle corresponding to the vector field of the FDE for a the
%   scalar time variable t and the state variable y.
%   The set of initial conditions gam is a vector of lenght K. 
%   eta(m) is equal to the vector of the collocation parameters. 
%   Constraint: 0 ≤ eta(1) ≤ · · · ≤eta(m).
%   r>=1 is grading exponent. For r=1 the mesh is uniform. For the optimal
%   setting of r, compare [1,3].
%   N is the number of mesh points.
%
%  IN OUTPUT
%   t: vector of the graded mesh: t(i+1)=b(i/N)^r, i=0,...,N.
%   y: vector of the approximate solution at the time steps t.
%   fval is the number of the function evaluation. 
%
%  REFERENCES
%
%  [1] A. Cardone, D. Conte, B. Paternoster, (2023). 
%      A MATLAB Code for Fractional Differential Equations Based on Two-Step 
%      Spline Collocation Methods. In: Cardone, A., Donatelli, M., Durastante, F., 
%      Garrappa, R., Mazza, M., Popolizio, M. (eds) Fractional Differential Equations. 
%      INDAM 2021. Springer INdAM Series, vol. 50, 121–146. Springer, Singapore.
%      https://doi.org/10.1007/978-981-19-7716-9_8
%
%  [2] A. Cardone, D. Conte, B. Paternoster, Stability of two-step spline collocation
%      methods for initial value problems for fractional differential equations,
%      Commun. Nonlinear Sci. Numer. Simul. 115, (2022), 106726.
%
%  [3] A. Cardone, D. Conte, B. Paternoster, Two-step collocation methods for 
%      fractional differential equations, Discrete Contin. Dyn. Syst. Ser. B 23 
%      (2018), no. 7, 2709–2725. 
%
%  [4] Cardone A., Conte D., Paternoster B. (2021) A MATLAB Implementation 
%      of Spline Collocation Methods for Fractional Differential Equations. 
%      Lect. Notes Comput. Sci., vol 12949, 387-401.
%
%%

t=b*((0:N)'/N).^r;
h=diff(t);
m=length(eta);
etatilde=linspace(0,1,2*m)';

% starting solution in [t(1),t(2)]
[Alagr_tilde]=matrix_Lagrange(etatilde);
[~,ytilde,ztilde,fval]=fcoll(f,t(2),gam,alpha,etatilde,r,1);
%disp('starting'), fval

H=[fliplr(vander(eta)),diag(eta.^m)*fliplr(vander(eta))];
Z=zeros(m,N); Z(:,1)=H*Alagr_tilde'*ztilde;
y=zeros(N+1,1); y(1)=gam(1); y(2)=ytilde(2);

Alagr=zeros(2*m,2*m,N);
options=optimset('TolFun',1e-14,'TolX',1e-14,'Display','off');
% options=optimset('TolFun',1e-6,'TolX',1e-6,'Display','off');

for j=2:N
    tj=t(j)+eta*h(j);
    Alagr(:,:,j)=tsmatrix_Lagrange(eta,h,j);
    A=tsmatrix_A(alpha,eta,Alagr(:,:,j));
    E=tsmatrix_E(alpha,eta,Alagr,t,h,j);
    Bj=tslag(Z,ztilde,A,E,Alagr,Alagr_tilde,t,h,eta,alpha,j);
    Qj=Q(tj,alpha,gam);
    iniz=Z(:,j-1);
    [Z(:,j),~,~,OUTPUT]=fsolve(@tssystem_F,iniz,options,f,tj,A,Bj,Qj,h,j,alpha);
    fval=fval+OUTPUT.funcCount;
    
    b=tsmatrix_A(alpha,1,Alagr(:,:,j));
    g=tsmatrix_E(alpha,1,Alagr,t,h,j);
    Wj=tslag(Z,ztilde,b,g,Alagr,Alagr_tilde,t,h,1,alpha,j);
    Qjp1=Q(t(j+1),alpha,gam);
    y(j+1)=h(j)^(alpha)*b(:,m+1:2*m)*Z(:,j)+Wj+Qjp1;
end

end