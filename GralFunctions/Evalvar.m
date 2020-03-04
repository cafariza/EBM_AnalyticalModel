function [Ux_exact,Ux1_exact,Ux2_exact,Ux3_exact,alpha_exact,alpha1_exact,alpha2_exact]= Evalvar(E,C,gama,wn,L,K,Lam,lw, H, Ki, Ux, alpha)
% This function computes modal transverse displacement U + macroscopic
% rotation plus their derivatives in an exact manner (using analytical
% formulations). The analytical expressions are given by means of symbolic
% functions.
wn1=wn;
%-------------------------------------------------%
%                  Derivatives U(x)
%-------------------------------------------------%
Ux1=diff(Ux,1); Ux2=diff(Ux,2); Ux3=diff(Ux,3); 
%-------------------------------------------------%
%                  Derivatives alpha(x)
%-------------------------------------------------%
alpha1=diff(alpha,1); alpha2=diff(alpha,2); 

x_t=0:lw:H;
X=x_t/L;
Ux_exact=zeros(length(X),length(E(1,:)));
Ux1_exact=zeros(length(X),length(E(1,:)));
Ux2_exact=zeros(length(X),length(E(1,:)));
Ux3_exact=zeros(length(X),length(E(1,:)));

alpha_exact=zeros(length(X),length(E(1,:))); 
alpha1_exact=zeros(length(X),length(E(1,:)));
alpha2_exact=zeros(length(X),length(E(1,:)));
for i=1:length(E(1,:))

a=E(4,i); b=E(3,i);
d=E(1,i); f=E(5,i);
g=E(2,i); h=E(6,i);
wn=wn1(i);

for k= 1: length(X)
%-------------------------------------------------%
%          Exact Solution of Ux and alpha + derivatives
%-------------------------------------------------%
Ux_num=double(subs(Ux(X(k))));
Ux_ult=double(subs(Ux(X(end)))); 

Ux_exact(k, i)= Ux_num/Ux_ult;

Ux1_num=double(subs(Ux1(X(k)))); Ux2_num=double(subs(Ux2(X(k)))); Ux3_num=double(subs(Ux3(X(k))));
Ux1_exact(k, i)=Ux1_num/Ux_ult;Ux2_exact(k, i)=Ux2_num/Ux_ult;Ux3_exact(k, i)=Ux3_num/Ux_ult;

alpha_num=double(subs(alpha(X(k)))); 
alpha_exact(k, i)= alpha_num/Ux_ult; 

alpha1_num=double(subs(alpha1(X(k)))); alpha2_num=double(subs(alpha2(X(k)))); 
alpha1_exact(k, i)= alpha1_num/Ux_ult; alpha2_exact(k, i)= alpha2_num/Ux_ult;

end
end
end