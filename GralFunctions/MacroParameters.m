%%Macroscopic Parameters
% HPDM IMPLEMENTATION 
%   NOTE: This code computes the frequencies and modal shapes of an
%   eqivalent beam model according to the macroscopic parameters of a unit
%   cell of a frame structure by computing the stiffnesses using MatLab
%   and Castem functions.

% By: Carolina FRANCO. carolina.franco@ifsttar.fr
% IFSTTAR-SV/SDOA
% DATE: AUGUST 2019
% MODIF: 22/11/2019


% Macroscopic Constants

        %(a) Linear masses
            Lam=LinearMass(af,aw,hf,hw,lf,lw,pe);
                        
          
        % (b)Bending stiffness
            %Inner Stiffness >>> EIw Ki
            %Global bending >>> EI Kgb
        % (c) Shear stiffness >>> K
        
            [Ka,Kw, Ki, Kgb]=stiff_K(af,aw,hf,hw,lf,lw,Ew,Ef,length(aw),length(af));
            Iz=Kgb/Ew; %Inertia global of the story : It wont be used in this code but in FEM_HM.m
            Sect1=Lam/pe; %Area of the story : It wont be used in this code but in FEM_HM.m
            Iy=sum(aw)*hf^3/12; %Inertia of the story : It wont be used in this code but in FEM_HM.m
            % --------------------------------------------------------------------------
            
            %% EXTRACTING K FROM CASTEM
%Warning : in this file a modified file to compute K has been replaced. Be
%sure to well correct the file recall in the function ExCastemKM
[K,cmdout2]= ExCastemKM(Ka,aw1,aw2,af1,af2,af,N,nw)
if nw==2
    K=Ka;
end 