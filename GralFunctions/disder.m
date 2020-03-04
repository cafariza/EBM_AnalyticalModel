function [Ux,Ux1,Ux2,Ux3,Ux4,Ux5,Ux6]= disder(L, E, lw, H,b1,b2,b3)

x_t=0:lw:H;
x=x_t/L;
X=x;
B1=-b1.^2;
B2=b2.^2;
B3=b3.^2;

Ux=zeros(length(x),length(E(1,:)));Ux1=zeros(length(x),length(E(1,:)));
Ux2=zeros(length(x),length(E(1,:)));Ux3=zeros(length(x),length(E(1,:)));
Ux4=zeros(length(x),length(E(1,:)));Ux5=zeros(length(x),length(E(1,:)));
Ux6=zeros(length(x),length(E(1,:)));
clear x
for i=1:length(E(1,:))
for g=1:length(X)
syms x;
sym Ux_sc(x);
sym gh(x);


% gh(x)=E(2,i)*cosh(b3(i)*x)+E(6,i)*sinh(b3(i)*x);
gh(x)=((E(3,i)*(B1(i)^2-B2(i)^2)*B3(i))/(B1(i)*b1(i)*(B2(i)^2-B3(i)^2)))/(cosh(pi()*b2(i)/2)*cosh(pi()*b3(i)/2)*B1(i)*(-B2(i)^2+B3(i)^2)+cos(pi()*b1(i)/2)*(cosh(pi()*b2(i)/2)*(-B1(i)^2+B2(i)^2)*B3(i)+cosh(pi()*b3(i)/2)*B2(i)*(B1(i)^2-B3(i)^2)))*(cosh(pi()*b2(i)/2)*B1(i)*(-B2(i)^2+B3(i)^2)*(cosh(b3(i)*x)*sin(pi()*b1(i)/2)*b1(i)+b3(i)*sinh(b3(i)*(pi()/2-x)))-cos(pi()*b1(i)/2)*(B2(i)*(B1(i)^2-B3(i)^2)*(cosh(b3(i)*x)*sinh(pi()*b2(i)/2)*b2(i)-b3(i)*sinh(b3(i)*(pi()/2-x)))+b3(i)*sinh(b3(i)*x)*cosh(pi()*b2(i)/2)*(-B1(i)^2+B2(i)^2)*B3(i)));
Ux_sc(x)=E(4,i)*cos(b1(i)*x)+E(3,i)*sin(b1(i)*x)+E(1,i)*cosh(b2(i)*x)+E(5,i)*sinh(b2(i)*x)+gh(x);

UxUlti=Ux_sc(X(end));
Ux_sc1=diff(Ux_sc)/UxUlti;
Ux_sc2=diff(Ux_sc1)/UxUlti;
Ux_sc3=diff(Ux_sc2)/UxUlti;
Ux_sc4=diff(Ux_sc3)/UxUlti;
Ux_sc5=diff(Ux_sc4)/UxUlti;
Ux_sc6=diff(Ux_sc5)/UxUlti;

Ux(g,i)=Ux_sc(X(g))/UxUlti;
Ux1(g,i)=Ux_sc1(X(g));
Ux2(g,i)=Ux_sc2(X(g));
Ux3(g,i)=Ux_sc3(X(g));
Ux4(g,i)=Ux_sc4(X(g));
Ux5(g,i)=Ux_sc5(X(g));
Ux6(g,i)=Ux_sc6(X(g));
end
plot(Ux,X); grid on;
hold on  
end

xlabel('Displacement U*','FontSize',12); ylabel('x*','FontSize',12);
legend(': Mode 1', ': Mode 2',  ': Mode 3', ': Mode 4');
title(' Modal Shapes ','FontSize',12);
end