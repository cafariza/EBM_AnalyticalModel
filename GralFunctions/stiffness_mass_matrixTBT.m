function [K_TBT, M_TBT, K_TBT1, M_TBT1]=stiffness_mass_matrixTBT()
% TIMOSHENKO BEAM
% 
 syms Kgb h Lam k x I
% 
syms alfa
alfa=Kgb/(k*h^2);
K_TBT =[  (12*Kgb)/(h^3*(12*alfa + 1)),                          (6*Kgb)/(h^2*(12*alfa + 1)), -(12*Kgb)/(h^3*(12*alfa + 1)),                          (6*Kgb)/(h^2*(12*alfa + 1));
  (6*Kgb)/(h^2*(12*alfa + 1)),   (4*Kgb*(36*alfa^2 + 6*alfa + 1))/(h*(12*alfa + 1)),  -(6*Kgb)/(h^2*(12*alfa + 1)), -(2*Kgb*(72*alfa^2 + 12*alfa - 1))/(h*(12*alfa + 1));
 -(12*Kgb)/(h^3*(12*alfa + 1)),                         -(6*Kgb)/(h^2*(12*alfa + 1)),  (12*Kgb)/(h^3*(12*alfa + 1)),                         -(6*Kgb)/(h^2*(12*alfa + 1));
  (6*Kgb)/(h^2*(12*alfa + 1)), -(2*Kgb*(72*alfa^2 + 12*alfa - 1))/(h*(12*alfa + 1)),  -(6*Kgb)/(h^2*(12*alfa + 1)),   (4*Kgb*(36*alfa^2 + 6*alfa + 1))/(h*(12*alfa + 1))];


M_TBT =[(Lam*h*(1680*alfa^2 + 294*alfa + 13))/(35*(12*alfa + 1)^2), (Lam*h^2*(1260*alfa^2 + 231*alfa + 11))/(210*(12*alfa + 1)^2),      (3*Lam*h*(560*alfa^2 + 84*alfa + 3))/(70*(12*alfa + 1)^2), -(Lam*h^2*(2520*alfa^2 + 378*alfa + 13))/(420*(12*alfa + 1)^2);
 (Lam*h^2*(1260*alfa^2 + 231*alfa + 11))/(210*(12*alfa + 1)^2),               (Lam*h^3)/120 + (Lam*h^3)/(840*(12*alfa + 1)^2),  (Lam*h^2*(2520*alfa^2 + 378*alfa + 13))/(420*(12*alfa + 1)^2),                (Lam*h^3)/(840*(12*alfa + 1)^2) - (Lam*h^3)/120;
  (3*Lam*h*(560*alfa^2 + 84*alfa + 3))/(70*(12*alfa + 1)^2), (Lam*h^2*(2520*alfa^2 + 378*alfa + 13))/(420*(12*alfa + 1)^2),     (Lam*h*(1680*alfa^2 + 294*alfa + 13))/(35*(12*alfa + 1)^2), -(Lam*h^2*(1260*alfa^2 + 231*alfa + 11))/(210*(12*alfa + 1)^2);
-(Lam*h^2*(2520*alfa^2 + 378*alfa + 13))/(420*(12*alfa + 1)^2),               (Lam*h^3)/(840*(12*alfa + 1)^2) - (Lam*h^3)/120, -(Lam*h^2*(1260*alfa^2 + 231*alfa + 11))/(210*(12*alfa + 1)^2),                (Lam*h^3)/120 + (Lam*h^3)/(840*(12*alfa + 1)^2)];
 

% M_TBT=[(Lam*h*(1/3*alfa^2 + 7/10*alfa + 13/35))/((12*alfa + 1)^2), (Lam*h^2*(1/24*alfa^2 + 11/120*alfa + 11/210))/((12*alfa + 1)^2),      (Lam*h*(1/6*alfa^2 + 3/10*alfa + 9/70))/((12*alfa + 1)^2), -(Lam*h^2*(1/24*alfa^2 +3/40*alfa + 13/420))/((12*alfa + 1)^2);
% (Lam*h^2*(1/24*alfa^2 + 11/120*alfa + 11/210))/((12*alfa + 1)^2), (Lam*h^3*(1/105+1/60*alfa+1/120*alfa^2))/((12*alfa + 1)^2),  (Lam*h^2*(1/24*alfa^2 + 3/40*alfa + 13/420))/((12*alfa + 1)^2),                -(Lam*h^3*(1/140+1/60*alfa+1/120*alfa^2))/((12*alfa + 1)^2);
%    (Lam*h*(1/6*alfa^2 + 3/10*alfa + 9/70))/((12*alfa + 1)^2), (Lam*h^2*(1/24*alfa^2 + 3/40*alfa + 13/420))/((12*alfa + 1)^2),     (Lam*h*(1/3*alfa^2 + 7/10*alfa + 13/35))/((12*alfa + 1)^2), (Lam*h^2*(1/24*alfa^2 + 11/120*alfa + 11/210))/((12*alfa + 1)^2);
%  -(Lam*h^2*(1/24*alfa^2 +3/40*alfa + 13/420))/((12*alfa + 1)^2),               -(Lam*h^3*(1/140+1/60*alfa+1/120*alfa^2))/((12*alfa + 1)^2), (Lam*h^2*(1/24*alfa^2 + 11/120*alfa + 11/210))/((12*alfa + 1)^2), (Lam*h^3*(1/105+1/60*alfa+1/120*alfa^2))/((12*alfa + 1)^2)];
%  

K_TBT1=[k/h, k/2,-k/h, k/2; k/2, Kgb/h+h*k/4, -k/2, -Kgb/h+h*k/4; -k/h, -k/2,k/h, -k/2; k/2, -Kgb/h+h*k/4, -k/2, Kgb/h+h*k/4];
M_TBT1=[Lam*h/3,0, Lam*h/6,0; 0, I*h/3, 0,I*h/6;Lam*h/6,0, Lam*h/3,0;0, I*h/6, 0,I*h/3];

end 
