% FRAME STRUCTURE: Simple Frame ---Going back to local scale
% HPDM IMPLEMENTATION 
%   NOTE: This code computes the frequencies and modal shapes of an
%   eqivalent beam model according to the macroscopic parameters of a unit
%   cell of a frame structure by computing the stiffnesses using MatLab
%   and Castem functions.

% By: Carolina FRANCO. carolina.franco@ifsttar.fr
% IFSTTAR-SV/SDOA
%% --------------------------------------------------------------------------
% 1 FRAME
% --------------------------------------------------------------------------

%% --------------------------------------------------------------------------
% INPUT DATA 1
% --------------------------------------------------------------------------
clc, clear, close all
tic

%Walls' Parameters
Ew=30000/(1000^2);                               % Young's modulus:MPa=MN/m^2
pe=2300/1000^3;                                  % concrete density (per unit volume):kg/m^3
poisson= 0.2;

hw=1000;                                         % wall width:mm
lw=3000;                                         % wall length:mm
aw1=115.66;
aw2=115.66;
aw=[aw1;aw2];                               % wall thickness:mm VECTOR
nw=length(aw);
%Floors' Parameters
Ef=30000/(1000^2);                               % Young's modulus:MPa
hf=1000;                                         % floor width:mm
af1=115.66;
af2=115.66;
af=[af1];                                    % floor thickness:mm VECTOR
lf=[3000];                                 % floor length:mm   VECTOR

% Macroscopic Constants

        %(a) Linear masses
            Lam=LinearMass(af,aw,hf,hw,lf,lw,pe)/(10^9);
                        
          
        % (b)Bending stiffness
            %Inner Stiffness >>> EIw Ki
            %Global bending >>> EI Kgb
        % (c) Shear stiffness >>> K
        
            [K, Ki, Kgb]=stiff_K(af,aw,hf,hw,lf,lw,Ew,Ef,length(aw),length(af));
            

% --------------------------------------------------------------------------
% INPUT DATA 2
% -------------------------------------------------------------------------- 

N=8;                   %Number of stories
d=4;                   %Number of modes
            
           H=lw*N;                  %Height of the structure mm

%      ffact=[1 1
%               3 6.25
%               5 17.36
%               7  34.03
%              ];             %Lower and Upper Limit of the frequency (Shear-Global Bending) [Hz]
%% EXTRACTING K FROM CASTEM
[K,cmdout2]= ExCastemKM(K,aw1,aw2,af1,af2,af,N,nw)


%% HPDM: NATURAL FREQUENCY COMPUTATIONS

omegar=zeros(d,1);
% --------------------------------------------------------------------------
% Global structural behavior
% --------------------------------------------------------------------------
x=zeros(d,1);
y=zeros(d,1);
C=zeros(d,1);
L=zeros(d,1);
epsilon=zeros(d,1);

for k=1:d
% Scale Factor
L(k)=2*N*lw/(pi*(2*k-1));
epsilon(k)=pi*(2*k-1)/(2*N);
% Macroscopic Constants
C(k)=Kgb/(K*(L(k))^2);
gamma=Ki/Kgb;
xE(k)=log(C(k))/log(epsilon(k));
yE(k)=log(gamma)/log(epsilon(k));
end

% --------------------------------------------------------------------------
% Frequency estimation
% --------------------------------------------------------------------------

omega=0:0.4:10000;
deter=zeros(length(omega),1); 
frange=cell(d,1,1);
% Scale Factor
L=2*N(1)*lw/(pi);
% Macroscopic Constants
C=Kgb/(K*(L)^2);
beta=C*(L/lw)^2;
gamma=Ki/Kgb;
for i=1:length(omega)
omegamsq1=Lam*omega(i)^2*(L)^2/K;

p=-((1+gamma)^2/(3*C^2*gamma^2)+omegamsq1/(C*gamma));
q=2*(1+gamma)^3/(27*C^3*gamma^3)-omegamsq1/(C^2*gamma)*(1-(1+gamma)/(3*gamma));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gamma)/(3*C*gamma);
B2=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gamma)/(3*C*gamma);
B3=2*sqrt(-p/3)*cos(eqn1)+(1+gamma)/(3*C*gamma);

b1=sqrt(-B1);
b2=sqrt(B2);
b3=sqrt(B3);
deter(i)=cos(pi*b1/2)*(B1^2-B2^2)*(B1^2-B3^2)*(2*B2*B3/(cosh(pi*b2/2)*cosh(pi*b3/2))+tanh(pi*b2/2)*tanh(pi*b3/2)*b2*b3*(B2+B3))+(B2^2-B1^2)*(B2^2-B3^2)*(2*B1*B3/(cosh(pi*b3/2))-sin(pi*b1/2)*tanh(pi*b3/2)*b1*b3*(B1+B3))+(B3^2-B1^2)*(B3^2-B2^2)*(2*B1*B2/(cosh(pi*b2/2))-sin(pi*b1/2)*tanh(pi*b2/2)*b1*b2*(B1+B2))-cos(pi*b1/2)*(B1^2*(B2^2-B3^2)^2+B3^2*(B2^2-B1^2)^2+B2^2*(B1^2-B3^2)^2);
end
n=0;
for i=1:length(deter)
    
    if deter(i)*deter(i+1)<0
        frange{i}=[omega(i)^2 omega(i+1)^2];
        n=n+1;
        
        if frange{i}==[0 0.1600]
         frange{i}=[];
         n=n-1;
        end
        
        switch n
            case d
                break
        end
    end
end
frange= frange(~cellfun('isempty',frange));

syms omegasq 
omegamsq=sqrt(Lam*omegasq*(L)^2/K);

p=-((1+gamma)^2/(3*C^2*gamma^2)+omegamsq^2/(C*gamma));
q=2*(1+gamma)^3/(27*C^3*gamma^3)-omegamsq^2/(C^2*gamma)*(1-(1+gamma)/(3*gamma));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gamma)/(3*C*gamma);
B2=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gamma)/(3*C*gamma);
B3=2*sqrt(-p/3)*cos(eqn1)+(1+gamma)/(3*C*gamma);

b1=sqrt(-B1);
b2=sqrt(B2);
b3=sqrt(B3);

eqn=cos(pi*b1/2)*(B1^2-B2^2)*(B1^2-B3^2)*(2*B2*B3/(cosh(pi*b2/2)*cosh(pi*b3/2))+tanh(pi*b2/2)*tanh(pi*b3/2)*b2*b3*(B2+B3))+(B2^2-B1^2)*(B2^2-B3^2)*(2*B1*B3/(cosh(pi*b3/2))-sin(pi*b1/2)*tanh(pi*b3/2)*b1*b3*(B1+B3))+(B3^2-B1^2)*(B3^2-B2^2)*(2*B1*B2/(cosh(pi*b2/2))-sin(pi*b1/2)*tanh(pi*b2/2)*b1*b2*(B1+B2))-cos(pi*b1/2)*(B1^2*(B2^2-B3^2)^2+B3^2*(B2^2-B1^2)^2+B2^2*(B1^2-B3^2)^2)==0;

for n=1:d
omegar(n)=vpasolve(eqn,omegasq,[frange{n}(1),frange{n}(2)],'random',true);
end

wn=sort(omegar);
Ndecimals = 5 ;
fround = 10.^Ndecimals; 
wn = round(fround*wn)/fround;
if wn(1)==0
wn(1) = [];
d=d-1;
end


    for s=1:d
    f_natural(s)=sqrt(wn(s))/(2*pi);
    f_ratio(s)=f_natural(s)/f_natural(1);
    end
    
f_natural    %Hz
f_ratio

% --------------------------------------------------------------------------
% Modal Displacement
% --------------------------------------------------------------------------
E=zeros(6,d);
b1=zeros(d,1);
b2=zeros(d,1);
b3=zeros(d,1);
omegam=zeros(d,1);
for i=1:d
omegamsqr=sqrt(Lam*wn(i)*(L)^2/K);
omegam(i)=omegamsqr^2;
p=-((1+gamma)^2/(3*C^2*gamma^2)+omegamsqr^2/(C*gamma));
q=2*(1+gamma)^3/(27*C^3*gamma^3)-omegamsqr^2/(C^2*gamma)*(1-(1+gamma)/(3*gamma));
eqn1=1/3*acos((q/2)*sqrt(-27/p^3));

B1=2*sqrt(-p/3)*cos(eqn1+2*pi/3)+(1+gamma)/(3*C*gamma);
B2=2*sqrt(-p/3)*cos(eqn1+4*pi/3)+(1+gamma)/(3*C*gamma);
B3=2*sqrt(-p/3)*cos(eqn1)+(1+gamma)/(3*C*gamma);

b1(i)=sqrt(-B1);
b2(i)=sqrt(B2);
b3(i)=sqrt(B3);

% cl1=t*(sinh(pi*b2/2)/b2-sin(pi*b1/2)/b1)+v*(sinh(pi*b3/2)/b3-sin(pi*b1/2)/b1)+s*(-cos(pi*b1/2)/b1+cosh(pi*b2/2)*(B2*(B1^2-B3^2))/(b1*B1*(B2^2-B3^2))+cosh(pi*b3/2)*(B3*(B2^2-B1^2))/(b1*B1*(B2^2-B3^2)))==0;
% cl2=t*(B2*cosh(pi*b2/2)-B1*cos(pi*b1/2))+v*(B3*cosh(pi*b3/2)-B1*cos(pi*b1/2))+s*(B1*sin(pi*b1/2)+sinh(pi*b2/2)*(b2*B2^2*(B1^2-B3^2))/(b1*B1*(B2^2-B3^2))+sinh(pi*b3/2)*(b3*B3^2*(B2^2-B1^2))/(b1*B1*(B2^2-B3^2)))==0;
% cl3=t*(cosh(pi*b2/2)*(B2+C*B1*B3)-cos(pi*b1/2)*(B1+C*B2*B3))+v*(cosh(pi*b3/2)*(B3+C*B1*B2)-cos(pi*b1/2)*(B1+C*B2*B3))+s*(sin(pi*b1/2)*(B1+C*B2*B3)+sinh(pi*b2/2)*(b2*B2*(B1^2-B3^2))/(b1*B1*(B2^2-B3^2))*(B2+C*B1*B3)-sinh(pi*b3/2)*(b3*B3*(B1^2-B2^2))/(b1*B1*(B2^2-B3^2))*(B3+C*B1*B2))==0;

t1=(sinh(pi*b2(i)/2)/b2(i)-sin(pi*b1(i)/2)/b1(i));
t2=(B2*cosh(pi*b2(i)/2)-B1*cos(pi*b1(i)/2));
t3=(cosh(pi*b2(i)/2)*(B2+C*B1*B3)-cos(pi*b1(i)/2)*(B1+C*B2*B3));

v1=(sinh(pi*b3(i)/2)/b3(i)-sin(pi*b1(i)/2)/b1(i));
v2=(B3*cosh(pi*b3(i)/2)-B1*cos(pi*b1(i)/2));
v3=(cosh(pi*b3(i)/2)*(B3+C*B1*B2)-cos(pi*b1(i)/2)*(B1+C*B2*B3));

s1=(-cos(pi*b1(i)/2)/b1(i)+cosh(pi*b2(i)/2)*(B2*(B1^2-B3^2))/(b1(i)*B1*(B2^2-B3^2))+cosh(pi*b3(i)/2)*(B3*(B2^2-B1^2))/(b1(i)*B1*(B2^2-B3^2)));
s2=(B1*sin(pi*b1(i)/2)+sinh(pi*b2(i)/2)*(b2(i)*B2^2*(B1^2-B3^2))/(b1(i)*B1*(B2^2-B3^2))+sinh(pi*b3(i)/2)*(b3(i)*B3^2*(B2^2-B1^2))/(b1(i)*B1*(B2^2-B3^2)));
s3=(sin(pi*b1(i)/2)*(B1+C*B2*B3)+sinh(pi*b2(i)/2)*(b2(i)*B2*(B1^2-B3^2))/(b1(i)*B1*(B2^2-B3^2))*(B2+C*B1*B3)-sinh(pi*b3(i)/2)*(b3(i)*B3*(B1^2-B2^2))/(b1(i)*B1*(B2^2-B3^2))*(B3+C*B1*B2));

A=[t1 v1 s1; t2 v2 s2; t3 v3 s3];
As=null(A);
E(1:3,i)=As(:,1);    % Matrix that contains t,v,s,r,u and w values for each mode of vibration
E(4,i)=-E(1,i)-E(2,i);
E(5,i)=-E(3,i)*b2(i)*B2*(B3^2-B1^2)/(b1(i)*B1*(B2^2-B3^2));
E(6,i)=E(3,i)*b3(i)*B3*(B2^2-B1^2)/(b1(i)*B1*(B2^2-B3^2));
end

C
xE(1)
yE(1)
gamma

%Displacement U*

% for i=1:length(E(1,:))
% Ux_sc=@(x) E(4,i)*cos(b1(i)*x)+E(3,i)*sin(b1(i)*x)+E(1,i)*cosh(b2(i)*x)+E(5,i)*sinh(b2(i)*x)+E(2,i)*cosh(b3(i)*x)+E(6,i)*sinh(b3(i)*x);
% end

% --------------------------------------------------------------------------
% Outputs: Plotting section
% --------------------------------------------------------------------------

figure 
plot(xE(1),yE(1),'or',xE(2),yE(2),'xb',xE(3),yE(3),'sg',xE(4),yE(4),'dm','LineWidth',1,'MarkerSize',6);
hold on
xlim([-3,3]); ylim([0,5]); gx=[-3 -1 0 1]; gy=[4 2 1 0];ix=[-1 -2 -3 -4]; iy=[0 1 2 3];
plot(xlim,  [1 1], '--k', [1 1], ylim,[-1 -1], ylim, gx, gy, ix, iy)           
hold off
legend(': Structure F1: mode 1', ': Structure F1: mode 2',': Structure F1: mode 3',': Structure F1: mode 4')
xlabel(' x (C= e^x) Bending/Shear','FontSize',12); ylabel('y (Gamma= e^y) Inner Bending/Global Bending','FontSize',12);
title(' Evolution of the structural behavior  ','FontSize',12);

figure 
loglog(N(1),f_natural(1),'sr',N(1),f_natural(2),'ok',N(1),f_natural(3),'dm','LineWidth',1,'MarkerSize',6); grid on;
hold on
xlabel('Cell number','FontSize',12); ylabel('Frequency [Hz]','FontSize',12);
legend(': Mode 1', ': Mode 2',  ': Mode 3',  ': Mode 4');
title(' Frequencies ','FontSize',12);

figure 
plot(omega,  deter, '-b'); grid on; 
xlabel('Omega [rad/s]','FontSize',12); ylabel('Amplitude','FontSize',12);

figure 
x_t=0:3000:H;
x=x_t/L;
Ux=zeros(length(x),length(E(1,:)));
Ux1=zeros(length(x),length(E(1,:)));
Ux2=zeros(length(x),length(E(1,:)));
Ux3=zeros(length(x),length(E(1,:)));
Ux4=zeros(length(x),length(E(1,:)));
Ux5=zeros(length(x),length(E(1,:)));
Ux6=zeros(length(x),length(E(1,:)));

Ux_sc=zeros(length(x),1);
Ux_sc1=zeros(length(x),1);
Ux_sc2=zeros(length(x),1);
Ux_sc3=zeros(length(x),1);
Ux_sc4=zeros(length(x),1);
Ux_sc5=zeros(length(x),1);
Ux_sc6=zeros(length(x),1);

for i=1:length(E(1,:))
for g=1:length(x)
UxUlti=E(4,i)*cos(b1(i)*x(end))+E(3,i)*sin(b1(i)*x(end))+E(1,i)*cosh(b2(i)*x(end))+E(5,i)*sinh(b2(i)*x(end))+E(2,i)*cosh(b3(i)*x(end))+E(6,i)*sinh(b3(i)*x(end));
Ux_sc(g)=(E(4,i)*cos(b1(i)*x(g))+E(3,i)*sin(b1(i)*x(g))+E(1,i)*cosh(b2(i)*x(g))+E(5,i)*sinh(b2(i)*x(g))+E(2,i)*cosh(b3(i)*x(g))+E(6,i)*sinh(b3(i)*x(g)))/UxUlti;
Ux_sc1(g)=(-E(4,i)*b1(i)*sin(b1(i)*x(g))+E(3,i)*b1(i)*cos(b1(i)*x(g))+E(1,i)*b2(i)*sinh(b2(i)*x(g))+E(5,i)*(i)*cosh(b2(i)*x(g))+E(2,i)*b3(i)*sinh(b3(i)*x(g))+E(6,i)*b3(i)*cosh(b3(i)*x(g)))/UxUlti;
Ux_sc2(g)=(-E(4,i)*b1(i)^2*cos(b1(i)*x(g))-E(3,i)*b1(i)^2*sin(b1(i)*x(g))+E(1,i)*b2(i)^2*cosh(b2(i)*x(g))+E(5,i)*b2(i)^2*sinh(b2(i)*x(g))+E(2,i)*b3(i)^2*cosh(b3(i)*x(g))+E(6,i)*b3(i)^2*sinh(b3(i)*x(g)))/UxUlti;
Ux_sc3(g)=(E(4,i)*b1(i)^3*sin(b1(i)*x(g))-E(3,i)*b1(i)^3*cos(b1(i)*x(g))+E(1,i)*b2(i)^3*sinh(b2(i)*x(g))+E(5,i)*b2(i)^3*cosh(b2(i)*x(g))+E(2,i)*b3(i)^3*sinh(b3(i)*x(g))+E(6,i)*b3(i)^3*cosh(b3(i)*x(g)))/UxUlti;
Ux_sc4(g)=(E(4,i)*b1(i)^4*cos(b1(i)*x(g))+E(3,i)*b1(i)^4*sin(b1(i)*x(g))+E(1,i)*b2(i)^4*cosh(b2(i)*x(g))+E(5,i)*b2(i)^4*sinh(b2(i)*x(g))+E(2,i)*b3(i)^4*cosh(b3(i)*x(g))+E(6,i)*b3(i)^4*sinh(b3(i)*x(g)))/UxUlti;
Ux_sc5(g)=(-E(4,i)*b1(i)^5*sin(b1(i)*x(g))+E(3,i)*b1(i)^5*cos(b1(i)*x(g))+E(1,i)*b2(i)^5*sinh(b2(i)*x(g))+E(5,i)*b2(i)^5*cosh(b2(i)*x(g))+E(2,i)*b3(i)^5*sinh(b3(i)*x(g))+E(6,i)*b3(i)^5*cosh(b3(i)*x(g)))/UxUlti;
Ux_sc6(g)=(-E(4,i)*b1(i)^6*cos(b1(i)*x(g))-E(3,i)*b1(i)^6*sin(b1(i)*x(g))+E(1,i)*b2(i)^6*cosh(b2(i)*x(g))+E(5,i)*b2(i)^6*sinh(b2(i)*x(g))+E(2,i)*b3(i)^6*cosh(b3(i)*x(g))+E(6,i)*b3(i)^6*sinh(b3(i)*x(g)))/UxUlti;
Ux(g,i)=Ux_sc(g);
Ux1(g,i)=Ux_sc1(g)/L;
Ux2(g,i)=Ux_sc2(g)/L^2;
Ux3(g,i)=Ux_sc3(g)/L^3;
Ux4(g,i)=Ux_sc4(g)/L^4;
Ux5(g,i)=Ux_sc5(g)/L^5;
Ux6(g,i)=Ux_sc6(g)/L^6;
end
plot(Ux_sc,x); grid on;
hold on  
end
xlabel('Displacement U*','FontSize',12); ylabel('x*','FontSize',12);
legend(': Mode 1', ': Mode 2',  ': Mode 3', ': Mode 4');
title(' Modal Shapes ','FontSize',12);

% --------------------------------------------------------------------------
% Modal Forces
% --------------------------------------------------------------------------
 alpha_sc=zeros(length(x),length(E(1,:)));
% T_sc=zeros(length(x),length(E(1,:)));
% M_sc=zeros(length(x),length(E(1,:)));
% Mi_sc=zeros(length(x),length(E(1,:)));
% 
% 
 for k=1:length(E(1,:))
     for g=1:length(x)
     alpha_sc(g,k)=(1+(Kgb/K^2)*Lam*wn(k))*Ux1(g,k)+(Kgb/K)*Ux3(g,k)-(Ki*Kgb/K^2)*Ux5(g,k);
     end
 end

% --------------------------------------------------------------------------
% Modal Forces
% --------------------------------------------------------------------------
[Theta_numuA, Theta_numuU] = ExCastemTheta(aw1,aw2,af1,af2,nw);
[T, M, thetan, thetaa]=forces_M(K,af,aw,hf,hw,lf,lw,Ew,Ux,Ux1,x,alpha_sc,pe,sqrt(wn),Theta_numuA,Theta_numuU,length(aw),length(af));
            
toc