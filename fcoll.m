function [t,y,Z,fval]=fcoll(f,b,gam,alpha,eta,r,N)
% one-step collocation method
% y solution at the mesh points
%
t=b*([0:N]'/N).^r;
h=diff(t);

Alagr=matrix_Lagrange(eta);
A=matrix_A(alpha,eta,Alagr);
bvect=vector_b(alpha,Alagr); % matrice A di dimensione s X m pagina 25 tesi federica laurino
options=optimset('TolFun',1e-14,'TolX',1e-14,'Display','off');
m=length(eta);
Z=zeros(m,N); y=zeros(N+1,1); y(1)=gam(1);
iniz=feval(f,t(1)+eta*h(1),gam(1)*ones(m,1));
fval=0;

for j=1:N
    tj=t(j)+eta*h(j);
    Bj=lag(Z,Alagr,t,h,eta,alpha,j);
    Qj=Q(tj,alpha,gam);
    [Z(:,j),~,~,OUTPUT]=fsolve(@system_F,iniz,options,f,tj,A,Bj,Qj,h,j,alpha);
    fval=fval+OUTPUT.funcCount;
    
    Wj=lag_y(Z,Alagr,t,h,alpha,j); 
    Qj=Q(t(j+1),alpha,gam);
    y(j+1)=h(j)^(alpha)*bvect'*Z(:,j)+Wj+Qj;
    
    iniz=Z(:,j);
end




end