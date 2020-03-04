function [N, T, M, M_in,M_in1, Tr,M_global,M_global1]=Internalforces(af,aw,hf,hw,lf,lw,Ew,Ux,Uy,Ux1,Ux2,Ux3, xd,alpha,alpha1,alpha2,pe,omega,nodal_rot,L, varargin)

%This function works for structures of Any number of frames. It computes
%the shear force and the bending moment of each structural member

%af         vector of floor thicknesses
%aw         vector of wall thicknesses
%hf         value of the width of the floor
%hw         value of the width of the wall
%lf         vector of length of floors
%lw         length of walls
%Ew         walls Young's modulus:MPa
%Ef         floor Young's modulus:MPa
%varargin   number of walls and number of floors

n=varargin{end-1};                               %Number of walls min 2
m=varargin{end};                                 %Number of floors min 1

[K,Kw, Ki, Kgb]=stiff_K(af,aw,hf,hw,lf,lw,Ew,Ew,n,m);
lcen=sum(lf)/2;
x=zeros(n,1);
d=zeros(n,1);
for i=1:m
    x(i+1,1)=lf(i)+x(i);
    d(i)=abs(lcen-x(i));
end
d(end)=d(1);
for i=1:n
Iw(i)=hw*aw(i)^3/12;                                   % moment of inertia:mm^4
Aw(i)=hw*aw(i);% cross section area of the wall:mm^2
end


l_T=zeros(n,length(omega));
l_L=zeros(n,length(omega));
l_Tas=zeros(n,length(omega));
l_Las=zeros(n,length(omega));
T={};
M={};

Ux(length(xd)+1,:)=zeros(1,length(omega));
for k=1:length(omega)
Uy {k}(length(xd)+1,:)=zeros(1,n);
end

Ux1(length(xd)+1,:)=zeros(1,length(omega));
%alpha=zeros(length(xd),length(omega));
alpha(length(xd)+1,:)=zeros(1,length(omega));
alpha1(length(xd)+1,:)=zeros(1,length(omega));
alpha2(length(xd)+1,:)=zeros(1,length(omega));
for j=1:length(omega)
for i=1:n
l_T(i,j)= (Ew*Iw(i)/((pe*Aw(i)*10^(-9))*(omega(j)^2)))^(1/4);
l_L(i,j)= (Ew/(pe*10^(-9)*(omega(j)^2)))^(1/2);
l_Tas(i,j)=lw/l_T(i,j);
l_Las(i,j)=lw/l_L(i,j);
for g=1:length(xd)

%For symmetric structures, the displacement at every node of a specific
%story is expected to be the same.
%Force in Newtons MN*10/10^6>>>N
Nw(g,i)=10^6*Ew*(Aw(i))/(l_L(i,j))*(Uy{j}(g,i)*cos(l_Las(i,j))-Uy{j}(g+1,i))/sin(l_Las(i,j));

Twtheta(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^3*(nodal_rot{j}(g,i)*sin(l_Tas(i,j))*sinh(l_Tas(i,j))-nodal_rot{j}(g+1,i)*(cos(l_Tas(i,j))-cosh(l_Tas(i,j))))*l_T(i,j)/(cos(l_Tas(i,j))*cosh(l_Tas(i,j))-1);
Tws(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^3*(Ux(g,j)*(cosh(l_Tas(i,j))*sin(l_Tas(i,j))+sinh(l_Tas(i,j))*cos(l_Tas(i,j)))-Ux(g+1,j)*(sin(l_Tas(i,j))+sinh(l_Tas(i,j))))/(cos(l_Tas(i,j))*cosh(l_Tas(i,j))-1);

        
Mwtheta(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^2*(nodal_rot{j}(g,i)*(cosh(l_Tas(i,j))*sin(l_Tas(i,j))-sinh(l_Tas(i,j))*cos(l_Tas(i,j)))-nodal_rot{j}(g+1,i)*(sin(l_Tas(i,j))-sinh(l_Tas(i,j))))*l_T(i,j)/(1-cos(l_Tas(i,j))*cosh(l_Tas(i,j)));
Mws(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^2*(Ux(g+1,j)*sin(l_Tas(i,j))*sinh(l_Tas(i,j))+Ux(g,j)*(cos(l_Tas(i,j))-cosh(l_Tas(i,j))))/(1-cos(l_Tas(i,j))*cosh(l_Tas(i,j)));
Mw(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^2*(Ux(g+1,j)*sin(l_Tas(i,j))*sinh(l_Tas(i,j))+Ux(g,j)*(cos(l_Tas(i,j))-cosh(l_Tas(i,j))))/(1-cos(l_Tas(i,j))*cosh(l_Tas(i,j)))+Mwtheta(g,i);
Tw(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^3*(Ux(g,j)*(cosh(l_Tas(i,j))*sin(l_Tas(i,j))+sinh(l_Tas(i,j))*cos(l_Tas(i,j)))-Ux(g+1,j)*(sin(l_Tas(i,j))+sinh(l_Tas(i,j))))/(cos(l_Tas(i,j))*cosh(l_Tas(i,j))-1)-Twtheta(g,i);
    if g==length(xd)
        Tw(g,i)=(Tw(g-1,i));
        Mw(g,i)=(Mw(g-1,i));
    end

Mi(g,i)=-Ki*Ux2(g,j)*10^6/(L^2);
Mi_1(g,i)=-Ki*Ux3(g,j)*10^6/(L^3);
Trs(g,i)=-K*(Ux1(g,j)/L-alpha(g,j))*10^6;
Mglobal(g,i)=-Kgb*alpha1(g,j)*10^6;
Mglobal1(g,i)=-Kgb*alpha2(g,j)*10^6;
end 
end
% Forces in the element
N{j}=[Nw];% Axial force
T{j}=[Tw];% Shear force
M{j}=[Mw];% Flexural moment

% From equilibrium equations
M_in{j}=Mi; % Inner bending moment
M_in1{j}=Mi_1;% First derivative inner bending moment
M_global{j}=Mglobal; % Global bending moment
M_global1{j}=Mglobal1; % First derivative global bending moment
Tr{j}=Trs; % Shear force
end
end