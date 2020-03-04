%% --------------------------------------------------------------------------
% Modal Displacement Coefficients
% --------------------------------------------------------------------------
[E, infsol]=Ecoeff(L,K,Lam, wn,C, gama,d);
% --------------------------------------------------------------------------
% Transverse Displacement, U and  Macroscopic Rotation, alpha
% -----------------------------------------------------------------------

x_t=0:lw:H;
x=x_t/L;

%Deformed shape is obtained from two alnertative solutions of the 6DOE
[Ux_a,Ux_b, alpha_a, alpha_b]= disdervar(); %Function that computes symbolic equations for U and alpha

%Solution a
[Ux,Ux1,Ux2,Ux3,alpha,alpha1,alpha2]= Evalvar(E,C,gama,wn,L,K,Lam,lw, H,Ki, Ux_b, alpha_b);