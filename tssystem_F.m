function F=tssystem_F(x,f,tj,A,Bj,Qj,h,j,alpha)
m=size(A,1);
F=x-feval(f,tj,h(j)^(alpha)*A(:,m+1:2*m)*x+Bj+Qj);

