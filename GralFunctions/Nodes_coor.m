function [XYZ]=Nodes_coor(lw,lf,N)
%% Create Model
%% Node coordinates 
% N=input('Number of cells: ');                  %Number of cells
%N=5:1:1000;   %Number of cells
Nf=N;
for j=1:length(N)
switch N(j)
case N
    div=lw;
otherwise
    div=lw/(N(j)/Nf);      
end
   
for i=1:2:(2*N(j)-1)
    XYZ(i,:)=[0 div*(i-1)/2];
end

for i=2:2:2*N(j)
    XYZ(i,:)= [lf XYZ(i-1,2)];
end 
end 
NUM_NOD=length(XYZ);
end