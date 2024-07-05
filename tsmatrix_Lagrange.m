function [Alagr]=tsmatrix_Lagrange(eta,h,j)

eta_h=[(eta'-1)*(h(j-1)/h(j)),eta']';
m_h=length(eta_h);
mu=1;
a=[1/(eta_h(mu)-eta_h(2)),-eta_h(2)/(eta_h(mu)-eta_h(2))];
for i=3:length(eta_h)
    b=[1/(eta_h(mu)-eta_h(i)),-eta_h(i)/(eta_h(mu)-eta_h(i))];
    a=conv(a,b);
end
Alagr(mu,:)=a;

for mu=2:m_h

    a=[1/(eta_h(mu)-eta_h(1)),-eta_h(1)/(eta_h(mu)-eta_h(1))];
    for i=2:mu-1
        b=[1/(eta_h(mu)-eta_h(i)),-eta_h(i)/(eta_h(mu)-eta_h(i))];
        a=conv(a,b);
    end
    for i=mu+1:m_h
        b=[1/(eta_h(mu)-eta_h(i)),-eta_h(i)/(eta_h(mu)-eta_h(i))];
        a=conv(a,b);
    end

    Alagr(mu,:)=a;
end

Alagr=fliplr(Alagr);
end
