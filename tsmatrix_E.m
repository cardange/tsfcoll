function E=tsmatrix_E(alpha,eta,Alagr,t,h,j)

s=length(eta);
m=size(Alagr,1)/2;
W=zeros(s,2*m);
E=zeros(s,2*m,j-1);
interv=0:2*m-1;

for lambda=2:j-1

    d=(t(j)-t(lambda)+eta*h(j))/h(lambda);
    for k=1:s
        W(k,:)= betainc(1/d(k),1+interv,alpha);
    end
    w=gamma(1+interv)./gamma(1+interv+alpha);
    E(:,:,lambda)=(((d.^(alpha+interv)).*w).*W)*Alagr(:,:,lambda)';

end

end