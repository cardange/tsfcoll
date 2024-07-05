
function lagterm=tslag(Z,ztilde,A,E,Alagr,Alagr_tilde,t,h,eta,alpha,j)

s=length(eta); m=size(Alagr,1)/2; N=length(t)-1;

hr=zeros(N,1);
for i=2:N
    hr(i)=h(i)/h(i-1);
end
hr=hr.^alpha;

lagterm=zeros(s,1);
if j==2
    B=hr(2)*A(:,1:m);
    lagterm=lagterm+h(1)^alpha*B*Z(:,1);
else
    B=hr(2)*E(:,1:m,2);
    lagterm=lagterm+h(1)^alpha*B*Z(:,1);
    for lambda=2:j-2
        B=E(:,m+1:2*m,lambda)+hr(lambda+1)*E(:,1:m,lambda+1);
        lagterm=lagterm+h(lambda)^alpha*B*Z(:,lambda);
    end
    B=hr(j)*A(:,1:m)+E(:,m+1:2*m,j-1);
    lagterm=lagterm+h(j-1)^alpha*B*Z(:,j-1);
end

interv=0:2*m-1;
d=(t(j)+eta*h(j))/h(1);
d2=(t(j)+eta*h(j));


W=zeros(s,2*m);
for k=1:s
    W(k,:)= betainc(1/d(k),1+interv,alpha);
end
Etilde=((d.^(alpha+interv)).*W)*Alagr_tilde';

Etilde=Etilde/gamma(alpha);
lagterm=lagterm+h(1)^alpha*Etilde*ztilde;

end