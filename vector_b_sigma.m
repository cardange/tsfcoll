function[b]=vector_b_sigma(alpha,Alagr,sigma)
% calcola vettore b della formula (20) dell'articolo Matlab implementation
% ma con sigma diverso da 1

m=size(Alagr,1);
interv=[0:m-1]';
w=gamma(interv+1)./gamma(interv+1+alpha);

b=Alagr*((sigma.^((alpha+interv))).*w);
end