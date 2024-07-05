clear
f=@(t,y) y.^2+1/gamma(1.5)*t.^0.5-t.^2;
b=1;
gam=[0];
alpha=1/2;
eta=[1/2,1]';
r=8;
N=b/2^-5;
[t,y]=tsfcoll(f,b,gam,alpha,eta,r,N);
err=abs(y(end)-b)

sol=@(t) t;
subplot(2,1,1)
plot(t,sol(t),'b',t,y,'*r','LineWidth',1.5)
legend('Exact solution','Numerical solution','Location','southeast')
xlabel('t'); ylabel('sol') 
set(gca,'Fontsize',12)

subplot(2,1,2)
semilogy(t,abs(y-sol(t)),'LineWidth',1.5)
%legend('Absolute error','Location','southeast')
grid
xlabel('t'); ylabel('abs err') 
set(gca,'Fontsize',12)


