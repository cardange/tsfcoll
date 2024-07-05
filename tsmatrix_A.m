function A=tsmatrix_A(alpha,eta,Alagr)

m=size(Alagr,1)/2;
interv=0:2*m-1;
w=gamma(interv+1)./gamma(interv+1+alpha);
A=((eta.^((alpha+interv))).*w)*Alagr';

end