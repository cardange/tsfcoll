function z=f7(t,y)

z=-y-y.^2+(1+erf(sqrt(t))).*exp(t)+exp(2*t);