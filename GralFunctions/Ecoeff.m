function [E, infsol]=Ecoeff(L,K,Lam, wn,C, gama,d)
%This function evaluates E vector (coefficients of the displacement
%solution of the 6 ODE of generic equation
%V[X_, \[Omega]_] =a[\[Omega]]*Cos[Subscript[\[Alpha], 1][\[Omega]]*X] + 
% b[\[Omega]]*Sin[Subscript[\[Alpha], 1][\[Omega]]*X] + 
% d[\[Omega]]*Cosh[Subscript[\[Alpha], 2][\[Omega]]*X] + 
% f[\[Omega]]*Sinh[Subscript[\[Alpha], 2][\[Omega]]*X] +
%+g*Cosh[Subscript[\[Alpha], 3]*X] + h*Sinh[Subscript[\[Alpha], 3]*X];

E=zeros(6,d);
B1=zeros(d,1);
B2=zeros(d,1);
B3=zeros(d,1);
b1=zeros(d,1);
b2=zeros(d,1);
b3=zeros(d,1);
a_co=zeros(1,d);
b_co=zeros(1,d);
d_co=zeros(1,d);
f_co=zeros(1,d);
g_co=zeros(1,d);
h_co=zeros(1,d);
for i=1:d
omegamsqr=sqrt(Lam*wn(i)*(L)^2/K);
p=-((1+gama)^2/(3*C^2*gama^2)+omegamsqr^2/(C*gama));
q=2*(1+gama)^3/(27*C^3*gama^3)-omegamsqr^2/(C^2*gama)*(1-(1+gama)/(3*gama));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1(i)=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gama)/(3*C*gama);
B2(i)=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gama)/(3*C*gama);
B3(i)=2*sqrt(-p/3)*cos(eqn1)+(1+gama)/(3*C*gama);
b1(i)=sqrt(-B1(i));
b2(i)=sqrt(B2(i));
b3(i)=sqrt(B3(i));

b_co(i)=1;
d_co(i)=((b_co(i)*(B1(i)^2-B3(i)^2)*B2(i))/(B1(i)*b1(i)*(B2(i)^2-B3(i)^2)))/(cosh(pi()*b2(i)/2)*cosh(pi()*b3(i)/2)*B1(i)*(-B2(i)^2+B3(i)^2)+cos(pi()*b1(i)/2)*(cosh(pi()*b2(i)/2)*(-B1(i)^2+B2(i)^2)*B3(i)+cosh(pi()*b3(i)/2)*B2(i)*(B1(i)^2-B3(i)^2)))*(cosh(pi()*b3(i)/2)*B1(i)*(B2(i)^2-B3(i)^2)*(sin(pi()*b1(i)/2)*b1(i)+sinh(pi()*b2(i)/2)*b2(i))+cos(pi()*b1(i)/2)*(B1(i)^2-B2(i)^2)*B3(i)*(sinh(pi()*b2(i)/2)*b2(i)-sinh(pi()*b3(i)/2)*b3(i)));
g_co(i)=((b_co(i)*(B1(i)^2-B2(i)^2)*B3(i))/(B1(i)*b1(i)*(B2(i)^2-B3(i)^2)))/(cosh(pi()*b2(i)/2)*cosh(pi()*b3(i)/2)*B1(i)*(-B2(i)^2+B3(i)^2)+cos(pi()*b1(i)/2)*(cosh(pi()*b2(i)/2)*(-B1(i)^2+B2(i)^2)*B3(i)+cosh(pi()*b3(i)/2)*B2(i)*(B1(i)^2-B3(i)^2)))*(cosh(pi()*b2(i)/2)*B1(i)*(-B2(i)^2+B3(i)^2)*(sin(pi()*b1(i)/2)*b1(i)+sinh(pi()*b3(i)/2)*b3(i))-cos(pi()*b1(i)/2)*(B1(i)^2-B3(i)^2)*B2(i)*(sinh(pi()*b2(i)/2)*b2(i)-sinh(pi()*b3(i)/2)*b3(i)));
a_co(i)=-d_co(i)-g_co(i);
f_co(i)=b_co(i)*b2(i)*B2(i)*(B1(i)^2-B3(i)^2)/(b1(i)*B1(i)*(B2(i)^2-B3(i)^2));
h_co(i)=b_co(i)*b3(i)*B3(i)*(B2(i)^2-B1(i)^2)/(b1(i)*B1(i)*(B2(i)^2-B3(i)^2));
end
E(1,1:d)=d_co;
E(2,1:d)=g_co;
E(3,1:d)=b_co;
E(4,1:d)=a_co;
E(5,1:d)=f_co;
E(6,1:d)=h_co;
infsol.b1=b1;
infsol.b2=b2;
infsol.b3=b3;
end