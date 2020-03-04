% FRAME STRUCTURE: MULTIPLE Frames ---Going back to local scale
% HPDM IMPLEMENTATION 
%   NOTE: This code computes the frequencies and modal shapes of the
%   eqivalent beam model (GENERIC) according to the macroscopic parameters of a unit
%   cell of a frame structure by computing the stiffnesses using MatLab
%   and Castem functions.

% By: Carolina FRANCO. carolina.franco@ifsttar.fr
% IFSTTAR-SV/SDOA
% DATE: AUGUST 2019
% MODIF: 22/11/2019
% MODIF: 03/04/2019
%% CLEANER of previous info
CleanStart
%% INPUT Section
Input
%% Start counting TIME
tic
%% HPDM: Lam, Ki, Kgb, K Estimation
tic
MacroParameters
toc
%% HPDM: NATURAL FREQUENCY COMPUTATIONS
tic 
FrequencyEstimation
toc
%% --------------------------------------------------------------------------
% MODAL DEFORMATION : Kinematic variables
% --------------------------------------------------------------------------
tic
KinematicVariables
toc
%% --------------------------------------------------------------------------
%NODAL ROTATION FROM HYBRID METHOD: ANALYTICAL + NUMERICAL 
% --------------------------------------------------------------------------
tic
[Theta_numuA, Theta_numuU] = ExCastemTheta(aw1,aw2,af1,af2,nw); %Unit Rotation from static analysis 
[UY_numuA, UY_numuU] =  ExCastemDisp(nw); %Unit Vertical displacement from static analysis 
toc
%% --------------------------------------------------------------------------
% Internal FORCES (ANALYTICAL RESULTS)
% --------------------------------------------------------------------------
tic
InForces
toc
%% End counting TIME
toc