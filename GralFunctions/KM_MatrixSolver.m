% This code generates the interpolation functions needed to compute the
% stiffness matrix and the mass matrix based on the Generic Formulation of
% the Equivalent Beam.

clear all
clc
%% Interpolation functions development
syms a b c d h
x=0;
for i=1:4
switch i
    case 1
        eqn1 = a+b*(x)+c*(x)^2+d*(x)^3  == 1;
eqn2 = b+2*c*(x)+3*d*(x)^2  == 0;
eqn3 = a+b*(x+h)+c*(x+h)^2+d*(x+h)^3  == 0;
eqn4 = b+2*c*(x+h)+3*d*(x+h)^2  == 0;
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [a, b, c, d])
X1 =  linsolve(A,B)
    case 2
        eqn1 = a+b*(x)+c*(x)^2+d*(x)^3  == 0;
eqn2 = b+2*c*(x)+3*d*(x)^2  == 1;
eqn3 = a+b*(x+h)+c*(x+h)^2+d*(x+h)^3  == 0;
eqn4 = b+2*c*(x+h)+3*d*(x+h)^2  == 0;
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [a, b, c, d])
X2=  linsolve(A,B)
    case 3
        eqn1 = a+b*(x)+c*(x)^2+d*(x)^3  == 0;
eqn2 = b+2*c*(x)+3*d*(x)^2  == 0;
eqn3 = a+b*(x+h)+c*(x+h)^2+d*(x+h)^3  == 1;
eqn4 = b+2*c*(x+h)+3*d*(x+h)^2  == 0;
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [a, b, c, d])
X3 =  linsolve(A,B)
    case 4
        eqn1 = a+b*(x)+c*(x)^2+d*(x)^3  == 0;
eqn2 = b+2*c*(x)+3*d*(x)^2  == 0;
eqn3 = a+b*(x+h)+c*(x+h)^2+d*(x+h)^3  == 0;
eqn4 = b+2*c*(x+h)+3*d*(x+h)^2  == 1;
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3, eqn4], [a, b, c, d])
X4 =  linsolve(A,B)
end
end 
%% Stiffness matrix coefficients' Computation (K)
% Computing first terms of the stiffness matrix
clear x
syms x
syms Ki Kgb k Lam
% Transverse Displacement: Polynomial functions Phi (1)
Phi11=1-3*(x/h)^2+2*x^3/h^3;
Phi12=x*(1-x/h)^2;
Phi13=3*(x/h)^2-2*(x/h)^3;
Phi14=x*((x/h)^2-x/h);
Phi1=[Phi11;Phi12;Phi13;Phi14];
Phi1t=Phi1.';
%Second derivatives of Phi (1)
Phi2d1=-6/h^2+12*x/h^3;
Phi2d2=-4/h+6*x/h^2;
Phi2d3=6/h^2-12*x/h^3;
Phi2d4=6*x/h^2-2/h;
Phi=[Phi2d1;Phi2d2;Phi2d3;Phi2d4]
Phit=Phi.'
%Stiffness Matrix of second derivative terms K11
for i=1:length(Phit)
    for j=1:length(Phi)
        K(i,j)=Ki*int(Phit(i)*Phi(j),[0,h]);
    end 
end

%First derivatives of Phi (1)
Phid1=-6*x/h^2+6*x^2/h^3;
Phid2=1-4*x/h+3*x^2/h^2;
Phid3=6*x/h^2-6*x^2/h^3;
Phid4=3*x^2/h^2-2*x/h;
Phid=[Phid1;Phid2;Phid3;Phid4]
Phidt=Phid.'
%Stiffness Matrix of first derivative terms K11
for i=1:length(Phidt)
    for j=1:length(Phid)
        Kd(i,j)=k*int(Phidt(i)*Phid(j),[0,h]);
    end 
end

K11=Kd+K% K11 = Sum of K and Kd (Inner bending and shear force due to transverse displacement)

% Alpha: Polynomial functions Phi (2)
Phi21=(1-x/h)*(1-2*x/h);
Phi22=4*x/h*(1-x/h);
Phi23=-x/h*(1-2*x/h);
Phi2=[Phi21;Phi22;Phi23]
Phi2t=Phi2.'
%First derivatives of Phi (2)
Phid21=diff(Phi21)
Phid22=diff(Phi22)
Phid23=diff(Phi23)
Phid2=[Phid21;Phid22;Phid23]
Phid2t=Phid2.'

%Stiffness Matrix of first derivative terms Phi(1) by function Phi(2) K12
for i=1:length(Phidt)
    for j=1:length(Phi2)
        K12(i,j)=-k*int(Phidt(i)*Phi2(j),[0,h]);
    end 
end

K21=K12.'

%Stiffness Matrix of first derivative terms of Phi(2) K22
for i=1:length(Phid2t)
    for j=1:length(Phid2)
        K22d(i,j)=Kgb*int(Phid2t(i)*Phid2(j),[0,h]);
    end 
end

%Stiffness Matrix of terms of Phi(2) K22
for i=1:length(Phi2t)
    for j=1:length(Phi2)
        K22(i,j)=k*int(Phi2t(i)*Phi2(j),[0,h]);
    end 
end

K22=K22d+K22;
KG(1:4,1:4)=K11;KG(1:4,5:7)=K12; KG(5:7,1:4)=K21; KG(5:7,5:7)=K22;
Ka=KG; Ka(:,6)=KG(:,7);Ka(:,7)=KG(:,6); 
Kaa=Ka;Kaa(6,:)=Kaa(7,:); Kaa(7,:)=Ka(6,:);
K_aa=Kaa(1:6,1:6);Kca=Kaa(7,1:6);
Kac=Kaa(1:6,7);Kcc=Kaa(7,7);Kred=K_aa-Kac*Kcc^(-1)*Kca;


%% Mass Matrix of terms of Phi(1) (M)
for i=1:length(Phi1t)
    for j=1:length(Phi1)
        M(i,j)=Lam*int(Phi1t(i)*Phi1(j),[0,h]);
    end 
end
%% SavinG DATA to export

save('Matrices','KG','Kred', 'M');
