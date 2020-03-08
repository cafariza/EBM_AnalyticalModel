% From equilibrium equations 

for j=1:length(omega)
for i=1:nw
for g=1:length(xd)
Mi(g,i)=-Ki*Ux2(g,j)*10^6/(L^2);
Mi_1(g,i)=-Ki*Ux3(g,j)*10^6/(L^3);
Trs(g,i)=-K*(Ux1(g,j)/L-alpha(g,j))*10^6;
Mglobal(g,i)=-Kgb*alpha1(g,j)*10^6;
Mglobal1(g,i)=-Kgb*alpha2(g,j)*10^6;
end
end

M_in{j}=Mi; % Inner bending moment
M_in1{j}=Mi_1;% First derivative inner bending moment
M_global{j}=Mglobal; % Global bending moment
M_global1{j}=Mglobal1; % First derivative global bending moment
Tr{j}=Trs; % Shear force
TotSF{j}=Tr{j}-M_in1{j};
end