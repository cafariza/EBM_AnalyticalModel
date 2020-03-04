% This program solves the <<Equivalent Beam>> problem based on the
% homogeneization method using Finite Element Method (FEM). 
CleanStart

%% INPUT DATA AND EXTERNAL DATA FROM MF_LocalScale.m and KM_MatrixSolver.m
run('KM_MatrixSolver')
%load('Matrices.mat')
clear all
run('MAIN')
%load('mainvar.mat')
close all
k=K;
load('Matrices.mat')
N_elem=input('How many elements do you want to use (Type 0 if you want to use the default number)?  ');
switch N_elem
    case 0
    N_elem=N; % Number of Finite Elements by default equal to number of stories
end
tic
%% FEM PROCEDURE FOR HM
Mele=M;
Kele=Kred;
ndg=3;
node=zeros(N_elem+1,2); % nodes
for i=1:N_elem+1
   node(i,1)=i; node(i,2)=H/N_elem*(i-1);
end
I=Iz;
[K_bc, M_bc]=amatrix(Kele,Mele, node, ndg, H,Ki,Kgb,k, Lam, I);

  %% EIGENVALUE PROBLEM : FUNDAMENTAL FRECUENCIES AND MODE SHAPES COMPUTATION
 
% --------------------------------------------------------------------------
% Calculate eigenvalues 
% --------------------------------------------------------------------------
%  
  ei=eig(K_bc,M_bc); % eigenvalues
% sorted natural angular frequencies [rad/s] 
  ef=sort(real(sqrt(ei))); 
% sorted natural frequencies [Hz]
  f_fem=ef/(2*pi);

 
% --------------------------------------------------------------------------
%  Mode shapes 
% --------------------------------------------------------------------------
%
[V,D]=eig(K_bc,M_bc); % eigenvalues and eigenvectors


toc
