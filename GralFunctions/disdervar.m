function [Ux_a,Ux_b, alpha_a, alpha_b]= disdervar()

syms x a b d f g h C wn L K Lam wn gama ;
sym Ux_a(x);
sym gh1(x);
sym Ux_sb(x);
sym gh2(x);

%-------------------------------------------------%
%                   Minor Variables
%-------------------------------------------------%

omegamsq=sqrt(Lam*wn*(L)^2/K);
omegam=omegamsq^2;
p=-((1+gama)^2/(3*C^2*gama^2)+omegamsq^2/(C*gama));
q=2*(1+gama)^3/(27*C^3*gama^3)-omegamsq^2/(C^2*gama)*(1-(1+gama)/(3*gama));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gama)/(3*C*gama);
B2=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gama)/(3*C*gama);
B3=2*sqrt(-p/3)*cos(eqn1)+(1+gama)/(3*C*gama);

b1=sqrt(-B1);
b2=sqrt(B2);
b3=sqrt(B3);
%-------------------------------------------------%
%                   Main Variables
%-------------------------------------------------%
gh1(x)=g*cosh(b3*x)+h*sinh(b3*x);
gh2(x)=((b*(B1^2-B2^2)*B3)/(B1*b1*(B2^2-B3^2)))/(cosh(pi*b2/2)*cosh(pi*b3/2)*B1*(-B2^2+B3^2)+cos(pi*b1/2)*(cosh(pi*b2/2)*(-B1^2+B2^2)*B3+cosh(pi*b3/2)*B2*(B1^2-B3^2)))*(cosh(pi*b2/2)*B1*(-B2^2+B3^2)*(cosh(b3*x)*sin(pi*b1/2)*b1+b3*sinh(b3*(pi/2-x)))-cos(pi*b1/2)*(B2*(B1^2-B3^2)*(cosh(b3*x)*sinh(pi*b2/2)*b2-b3*sinh(b3*(pi/2-x)))+b3*sinh(b3*x)*cosh(pi*b2/2)*(-B1^2+B2^2)*B3));
Ux_a(x)=a*cos(b1*x)+b*sin(b1*x)+d*cosh(b2*x)+f*sinh(b2*x)+gh1(x);
Ux_b(x)=a*cos(b1*x)+b*sin(b1*x)+d*cosh(b2*x)+f*sinh(b2*x)+gh2(x);

Ux_a1=diff(Ux_a); Ux_b1=diff(Ux_b);
Ux_a2=diff(Ux_a1); Ux_b2=diff(Ux_b1);
Ux_a3=diff(Ux_a2); Ux_b3=diff(Ux_b2);
Ux_a4=diff(Ux_a3); Ux_b4=diff(Ux_b3);
Ux_a5=diff(Ux_a4); Ux_b5=diff(Ux_b4);

alpha_a=((1+C*omegam)*Ux_a1+C*Ux_a3-C^2*gama*Ux_a5)/(L);
alpha_b=((1+C*omegam)*Ux_b1+C*Ux_b3-C^2*gama*Ux_b5)/(L);


end