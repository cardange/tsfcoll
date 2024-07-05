function Wj=lag_y_sigma(Z,Alagr,t,h,alpha,j,sigma)
% calcola matrice g della formula (21) dell'articolo con sigma qualunque
% in lag_y invece si considera sigma=1
Wj=0;
m=size(Alagr,1);
interv=[0:m-1]';
w=gamma(interv+1)./gamma(interv+1+alpha);
for lambda=1:j-1
    d=(t(j)+sigma*h(j)-t(lambda))/h(lambda);
    W=betainc(1/d,1+interv,alpha);
    g=Alagr*(((d.^((alpha+interv))).*w).*W);
    Wj=Wj+(h(lambda)^alpha)*g'*Z(:,lambda);
end
end