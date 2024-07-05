function z=f5(t,y)
global alpha5
alpha=alpha5;
z=-(1+t.^2).*y.^2+...
    (t.^(1-alpha))/((1-alpha)*gamma(1-alpha))+(1+t.^2).*(1+t).^2;