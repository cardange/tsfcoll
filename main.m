%%
%  Test problems presented in
%  [1] A. Cardone, D. Conte, B. Paternoster, (2023).
%      A MATLAB Code for Fractional Differential Equations Based on Two-Step
%      Spline Collocation Methods. In: Cardone, A., Donatelli, M., Durastante, F.,
%      Garrappa, R., Mazza, M., Popolizio, M. (eds) Fractional Differential Equations.
%      INDAM 2021. Springer INdAM Series, vol. 50, 121–146. Springer, Singapore.
%      https://doi.org/10.1007/978-981-19-7716-9_8
%%

clear
problem = 2;
eta=[0.5,1]';
%eta=[0;0.5];
if problem ==1
    % test n.1 [1]
    b=1;
    alpha=1/2;
    gam=0;
    r=8;
    f=@f1;
    sol=@sol1;
elseif problem ==2
    % test problem n.2 [1]
    b=1;
    alpha=0.5;
    f=@(t,y) 40320/gamma(9-alpha)*t.^(8-alpha)-3*gamma(5+alpha/2)/gamma(5-alpha/2)*t.^(4-alpha/2)+...
        9/4*gamma(alpha+1)+(3/2*t.^(alpha/2)-t.^4).^3-y.^(3/2);
    gam=[0];
    v=0.5;
    r=4/v;
    sol=@(t) t.^8-3*t.^(4+alpha/2)+9/4*t.^(alpha);
elseif  problem==3
    % test problem n.3 [1]
    b=8;
    global alpha3 lambda3 rho3 y03
    alpha3=0.3; lambda3=-3; rho3=0.8;
    alpha=alpha3;
    gam=[2];
    y03=gam;
    %r6=1.4; %c=[0 0.5]
    r=26.7;
    f=@f3;
    sol=@sol3;
elseif problem ==4
    % test problem n.4 [1]
    b=1;
    alpha=2.5;
    gam=[1;1;1];
    r=8;
    f=@f4;
    sol=@(t) exp(t);
elseif problem ==5
    % test problem n.5 [1]
    b=1;
    global alpha5
    alpha5=0.4;
    alpha=alpha5;
    gam=[1];
    y05=gam;
    %r6=1.4; %c=[0 0.5]
    r=4/(1-alpha);
    f=@f5;
    sol=@(t) 1+t;
end
N=b/(2^-5);
[t,y,fval]=tsfcoll(f,b,gam,alpha,eta,r,N);
err=abs(y(end)-sol(b));
rel_err=abs(y(end)-sol(b))/sol(b);
T = table(err, rel_err, fval);
format short g
disp(T)
format short

subplot(2,1,1)
t1=linspace(0,b,1000);
plot(t1,sol(t1),'b',t,y,'*r','LineWidth',1.5)
xlabel('t');ylabel('sol');
legend('Exact solution','Numerical solution','Location','southeast')
set(gca,'Fontsize',12)

subplot(2,1,2)
semilogy(t,abs(y-sol(t)),'LineWidth',1.5)
grid
xlabel('t');ylabel('abs error');
%legend('Absolute error','Location','northeast')
set(gca,'Fontsize',12)

