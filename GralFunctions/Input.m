% FRAME STRUCTURE: MULTIPLE Frames ---Going back to local scale
% HPDM IMPLEMENTATION 
%   NOTE: This code computes the frequencies and modal shapes of an
%   eqivalent beam model according to the macroscopic parameters of a unit
%   cell of a frame structure by computing the stiffnesses using MatLab
%   and Castem functions.

% By: Carolina FRANCO. carolina.franco@ifsttar.fr
% IFSTTAR-SV/SDOA
% DATE: AUGUST 2019
% MODIF: 22/11/2019
%% --------------------------------------------------------------------------
% N FRAMES
% --------------------------------------------------------------------------

 %--------------------------------------------------------------------------
% INPUT DATA 
% --------------------------------------------------------------------------



%Walls' Parameters
Ew=30000/(1000^2);                               % Young's modulus:MPa=MN/m^2
pe=2.300/10^15;                                  % concrete density (per unit volume):kg/m^3
poisson= 0.2;
%Ew=Ew/(1-poisson^2);

hw=1000;%<<<<<<<<<<<<<<                          % wall width:mm
lw=3000;%<<<<<<<<<<<<<<                          % wall length:mm
%Floors' Parameters
Ef=30000/(1000^2);                               % Young's modulus:MPa
hf=1000;%<<<<<<<<<<<<<<                          % floor width:mm

nw=input('How many walls?  ');
N=input('Number of stories: ');                   %Number of stories
d=4;                                                %Number of modes
aw=zeros(nw,1);                                  % wall thickness:mm VECTOR

for i=1:nw
    fit=num2str(i);
    confit=strcat('Thickness of wall ',fit,':');
    aw(i)=input(confit);
end
tf=input('Thickness of the slab: ');

%% Assignments
aw1=aw(1);%<<<<<<<<<<<<<<
aw2=aw(2);%<<<<<<<<<<<<<<
   
af1=tf;%<<<<<<<<<<<<<<
af2=tf;%<<<<<<<<<<<<<<

n=nw-1;                                          % Number of frames
af=zeros(n,1);
lf=zeros(n,1);
for i=1:n
af(i)=tf;                                         % floor thickness:mm VECTOR
lf(i)=lw;                                        % floor length:mm   VECTOR
end
H=lw*N;                                  %Height of the structure mm



