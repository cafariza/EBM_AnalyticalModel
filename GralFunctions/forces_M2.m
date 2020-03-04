function [T, M, M_in,M_in1,M_global, Tr]=forces_M2(af,aw,hf,hw,lf,lw,Ew,Ux,Ux1,Ux2,Ux3, xd,pe,omega,L,thetanum,alphanum,alphanum1, varargin)
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

nodal_rot=thetanum;
n=varargin{end-1};                               %Number of walls min 2
m=varargin{end};                                 %Number of floors min 1
[K,Kw, Ki, Kgb]=stiff_K(af,aw,hf,hw,lf,lw,Ew,Ew,n,m);



for i=1:n
Iw(i)=hw*aw(i)^3/12;                                   % moment of inertia:mm^4
Aw(i)=hw*aw(i);% cross section area of the wall:mm^2
end


  
l_T=zeros(n,length(omega));
l_L=zeros(n,length(omega));
l_Tas=zeros(n,length(omega));
T={};
M={};
for k=1:length(omega)
Ux {k}(length(xd)+1,:)=zeros(1,n);
Ux1 {k}(length(xd),:)=zeros(1,n);
Ux2 {k}(length(xd)-1,:)=zeros(1,n);
Ux2 {k}(length(xd),:)=zeros(1,n);
Ux3 {k}(length(xd)-2,:)=zeros(1,n);
Ux3 {k}(length(xd)-1,:)=zeros(1,n);
Ux3 {k}(length(xd),:)=zeros(1,n);
alphanum1 {k}(length(xd),:)=zeros(1,n);
alphanum1 {k}(length(xd)+1,:)=zeros(1,n);
nodal_rot {k}(length(xd)+1,:)=zeros(1,n);
end

for j=1:length(omega)
for i=1:n
l_T(i,j)= (Ew*Iw(i)/((pe*Aw(i)*10^(-9))*(omega(j)^2)))^(1/4);
l_L(i,j)= (Ew/(pe*(omega(j)^2)))^(1/4);
l_Tas(i,j)=lw/l_T(i,j);
l_Las(i,j)=lw/l_L(i,j);
for g=1:length(xd)

%For symmetric structures, the displacement at every node of a specific
%story is expected to be the same.
%Force in Newtons MN*10/10^6>>>N
%Nw=10^6*Ew*Aw(i)*10^-9)/(l_L(i,j))*(Uy(g,j)*cos(l_Las(i,j))-Uy(g+1,j))/sin(l_Las(i,j));

Twtheta(g,i)=(nodal_rot{j}(g,i)*sin(l_Tas(i,j))*sinh(l_Tas(i,j))-nodal_rot{j}(g+1,i)*(cos(l_Tas(i,j))-cosh(l_Tas(i,j))))*l_T(i,j);
Mwtheta(g,i)=(nodal_rot{j}(g,i)*(cosh(l_Tas(i,j))*sin(l_Tas(i,j))-sinh(l_Tas(i,j))*cos(l_Tas(i,j)))-nodal_rot{j}(g+1,i)*(sin(l_Tas(i,j))-sinh(l_Tas(i,j))))*l_T(i,j);
Tw(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^3*(Ux{j}(g,i)*(cosh(l_Tas(i,j))*sin(l_Tas(i,j))+sinh(l_Tas(i,j))*cos(l_Tas(i,j)))-Ux{j}(g+1,i)*(sin(l_Tas(i,j))+sinh(l_Tas(i,j)))-Twtheta(g,i))/(cos(l_Tas(i,j))*cosh(l_Tas(i,j))-1);
Mw(g,i)=10^6*Ew*Iw(i)/(l_T(i,j))^2*(Ux{j}(g+1,i)*sin(l_Tas(i,j))*sinh(l_Tas(i,j))+Ux{j}(g,i)*(cos(l_Tas(i,j))-cosh(l_Tas(i,j)))+Mwtheta(g,i))/(1-cos(l_Tas(i,j))*cosh(l_Tas(i,j)));
  if g==length(xd)
        Tw(g,i)=(Tw(g-1,i));
        Mw(g,i)=(Mw(g-1,i));
    end
Trs(g,i)=-K*(Ux1{j}(g,i)/(L^1)-alphanum{j}(g,i));
Mi(g,i)=-Ki*Ux2{j}(g,i)*10^6/(L^2);
Mi_1(g,i)=-Ki*Ux3{j}(g,i)*10^6/(L^3);
Mglobal(g,i)=-Kgb*alphanum1{j}(g,i)*10^6;
end 
end
T{j}=Tw;
M{j}=Mw;
M_in{j}=Mi;
M_in1{j}=Mi_1;
M_global{j}=Mglobal;
Tr{j}=Trs;
end