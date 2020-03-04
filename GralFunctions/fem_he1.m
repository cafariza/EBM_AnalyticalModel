% This program solves the <<Equivalent Beam>> problem based on the
% homogeneization method using Finite Element Method (FEM). 
clc
clear all
close all

%% INPUT DATA AND EXTERNAL DATA FROM MF_LocalScale.m and KM_MatrixSolver.m
%run('KM_MatrixSolver')
load('Matrices.mat')
clear all
%run('MF_LocalScale')
load('mainvar.mat')
close all
k=K;
load('Matrices.mat')
%% --------------------------------------------------------------------------
% Compare the first six modal testing and F.E.M values of frequencies with 
% analytical and numerical (CASTEM) values.
% --------------------------------------------------------------------------
% first three frequencies of Castem 
f_cast= [fscanf(fopen('FREQ1.inp','r'),'%f') ...
       fscanf(fopen('FREQ2.inp','r'),'%f') ...
       fscanf(fopen('FREQ3.inp','r'),'%f') ...
       ];
   fclose('all');
% first three frequencies of analytical results
f_an= f_natural(1,1:3)';
ModeShapes_cast=[Unum{1,1}(:,1),Unum{1,2}(:,1),Unum{1,3}(:,1)];
ModeShapes_an=Ux(:,1:3);

Ne_max=100;
r_ex=zeros(1,3);
r=zeros(Ne_max,3);
error_fc=zeros(1,3);
x_f=H/L;x_0=0;
dx1=lw/L;           %Element length : cell size of the HM beam and castem stucture
x_hm=[x_0:dx1:x_f];   %Position vector of the HM beam and castem stuctur
savemodes_he={};

C = {'c','m','r','g','y',[.5 .6 .7],[.8 .2 .6]}; % Cell array of colors.

for i=1:1:Ne_max
     N_elem=i;
     
% FEM PROCEDURE FOR HM
Mele=M;
Kele=Kred;
ndg=3;
node=zeros(N_elem+1,2); % nodes
for s=1:N_elem+1
   node(s,1)=s; node(s,2)=H/N_elem*(s-1);
end
I=Iz;
[K_bc, M_bc]=amatrix(Kele,Mele, node, ndg, H,Ki,Kgb,k, Lam, I);
%% EIGENVALUE PROBLEM : FUNDAMENTAL FRECUENCIES AND MODE SHAPES COMPUTATION
% --------------------------------------------------------------------------
% Calculate eigenvalues 
% --------------------------------------------------------------------------
  ei=eig(K_bc,M_bc); % eigenvalues
% sorted natural angular frequencies [rad/s] 
  ef=sort(real(sqrt(ei))); 
% sorted natural frequencies [Hz]
  f_fem=ef/(2*pi);
  for j=1:3
    error_fc(1,j)=100*(f_fem(j)-f_cast(j))/f_cast(j);
end
r(i,:)=error_fc;
[V,D]=eig(K_bc,M_bc); % eigenvalues and eigenvectors
[d,ind] = sort(diag(D));
Ds = D(ind,ind);
Vs= zeros(size(V)+3);
Vs(4:end,1:end-3) = V(:,ind);

eigval = D(1,1);
eigvec = V(:,1);
error=K_bc*eigvec - eigval*M_bc*eigvec;

v1 = Vs(1:(ndg):end,1:end) ;    % Collecting only displacement degree's of Freedom
%
Vec = zeros(size(v1)) ;
for l = 1:size(v1,1)
    for p=1:size(v1,2)
        v(l,p)=v1(l,p)-v1(1,p);
    end
end

for p = 1:size(v1,2)
   ma= abs(max(v(:,p)));
   mi= abs(min(v(:,p)));
    Vec(:,p) = v(:,p)./(max(ma,mi)) ;
end

%Location of analysed nodes
dx=H/(L*N_elem);    %Element length of the Discretized beam according to N_elem
x_efem=x_0:dx:x_f;       %Position vector of the Discretized beam according to N_elem

ModeShapes_fem=Vec(:,1:3); 

savemodes_he{1,i}=ModeShapes_fem;
end

for i=1:3
    r_ex(i)=100*(abs(f_an(i)-f_cast(i)))/f_cast(i);
end
figure(1); 
plot(1:1:3,r_ex,'r*-',1:1:3,r,'bo--','LineWidth', 2)
title('Relative errors: CASTEM as reference[%] ');
xlabel('Mode');
ylabel('[%]');
axis([1,3,0,10])
xbounds = xlim;
set(gca,'XTick',xbounds(1):xbounds(2));
for i=1:1:Ne_max
    if i>10 && i <100
        text(2.5,r(i,end)*0.96,' ');
    else
str = sprintf('   Ne = %d',i);
text(2.5,r(i,end)*0.96,str);
hold on
    end
end

