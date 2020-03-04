function [K,Kw, Ki, Kgb]=stiff_K(af,aw,hf,hw,lf,lw,Ew,Ef,varargin)

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
Aw(i)=hw*aw(i);                                        % cross section area of the wall:mm^2
I(i)=Aw(i)*(d(i))^2;
end
for i=1:m
If(i)=hf*af(i)^3/12;                                   % moment of inertia:mm^4
Af(i)=hf*af(i);                                        % cross section area of the floor:mm^2
end
% (b)Bending stiffness
%Inner Stiffness >>> EIw
Ki=sum(Ew.*Iw);

%Global bending >>> EI
Kgb=sum(Ew.*I);

  kw=zeros(n,1);
  for i=1:n 
    kw(i)=12*Ew*Iw(i)/(lw^3);
  end
  Kw=sum(kw)*lw; 
switch m
    case 1
% (c) Shear stiffness

  
  kf=12*Ef*If/(lf^3);
  Kf=kf*lf^2/lw;
  K=Kf*Kw/(Kw+Kf);
disp('K is analytical')
    case 2
% (c) Shear stiffness
    
   kf1=12*Ef*If(1)/(lf(1)*lw);
   kf2=12*Ef*If(2)/(lf(2)*lw);
   %K=1/((1/(kf1+kf2))+1/(Kw));         % Stafford Smith et all formula
   K=kf1*(kf1*(Kw)+6*((kw(1)*lw+kw(3)*lw)*kw(2)*lw))/((kf1^2)+2*kf1*(Kw)+3*((kw(1)*lw+kw(3)*lw)*kw(2)*lw)); %HPDM formula
   disp('K is analytical')
    otherwise
      K=0;
      disp('K was extracted from CASTEM')
end


end