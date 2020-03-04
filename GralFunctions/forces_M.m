function [N, T, M, M_in,M_in1, Tr,M_global,M_global1, thetaa, thetan]=forces_M(af,aw,hf,hw,lf,lw,Ew,Ux,Uy,Ux1,Ux2,Ux3, xd,alpha,alpha1,alpha2,pe,omega,Theta_numuA,Theta_numuU,L, varargin)

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
I(i)=Aw(i)*(d(i))^2;
end
for i=1:m
If(i)=hf*af(i)^3/12;                                   % moment of inertia:mm^4
Af(i)=hf*af(i);                                        % cross section area of the floor:mm^2
end

  
  
l_T=zeros(n,length(omega));
l_L=zeros(n,length(omega));
l_Tas=zeros(n,length(omega));
l_Las=zeros(n,length(omega));
T={};
M={};
thetaa={};
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
switch m
    case 1
    %The following equations are applicable only to single frame structures.
    %Here thetaa{j} has the same value for both structural members at each story.
    kf=12*Ew*If(1)/(lf(1)^3); %WARNING It's computed only for a SINGLE FRAME HERE
    Kf=kf*lf(1)^2/lw;
  
    thetaa{j}(g,i)=K*(alpha(g,j)/Kw+(Ux1(g,j)/(Kf*L)))*-1;
    thetaa{j}(g+1,i)=K*(alpha(g+1,j)/Kw+(Ux1(g+1,j)/(Kf*L)))*-1;

    thetan{j}(g,i)=(alpha(g,j)*Theta_numuA(i))+(Ux1(g,j)*Theta_numuU(i)/(1/(lw/L)));
    thetan{j}(g+1,i)=(alpha(g+1,j)*Theta_numuA(i))+(Ux1(g+1,j)*Theta_numuU(i)/(1/(lw/L)));

    nodal_rot{j}(g,i)=thetaa{j}(g,i);
    nodal_rot{j}(g+1,i)=thetaa{j}(g+1,i);

    if j==length(omega)
    disp('Theta is computed by both ana and num. SF');
    end
    case 2
    %For double frame buildings 
    thetan{j}(g,i)=(alpha(g,j)*Theta_numuA(i))+(Ux1(g,j)*Theta_numuU(i)/(1/(lw/L)));
    thetan{j}(g+1,i)=(alpha(g+1,j)*Theta_numuA(i))+(Ux1(g+1,j)*Theta_numuU(i)/(1/(lw/L)));

    nodal_rot{j}(g,i)=thetan{j}(g,i);
    nodal_rot{j}(g+1,i)=thetan{j}(g+1,i);
    if j==length(omega)
    disp('Theta is computed numerically. DF');
    end
    otherwise 
    %Multiple Frame
    %Theta must be computed from the superposition of U_a*Theta_num(unitU) and
    %Alpha_a*Theta_num(unitAlpha) 
    thetan{j}(g,i)=(alpha(g,j)*Theta_numuA(i))+(Ux1(g,j)*Theta_numuU(i)/(1/(lw/L)));
    thetan{j}(g+1,i)=(alpha(g+1,j)*Theta_numuA(i))+(Ux1(g+1,j)*Theta_numuU(i)/(1/(lw/L)));

    nodal_rot{j}(g,i)=thetan{j}(g,i);
    nodal_rot{j}(g+1,i)=thetan{j}(g+1,i);
    if j==length(omega)
    disp('Theta is computed numerically');
    end
end 
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
N{j}=[Nw];
T{j}=[Tw];
M{j}=[Mw];


M_in{j}=Mi;
M_in1{j}=Mi_1;
M_global{j}=Mglobal;
M_global1{j}=Mglobal1;
Tr{j}=Trs;
end
end