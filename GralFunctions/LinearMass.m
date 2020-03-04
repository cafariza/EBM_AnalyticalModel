function [Lam]=LinearMass(af,aw,hf,hw,lf,lw,pe)

%af         vector of floor thicknesses
%aw         vector of wall thicknesses
%hf         value of the width of the floor
%hw         value of the width of the wall
%lf         vector of length of floors
%lw         length of walls
%pe         Material density

Lamw=pe*sum(hw.*aw)*lw;                                          % Mass of the walls
Lamf=pe*sum(hf.*af.*lf);                                         % Mass of the floors
Lam=(Lamw+Lamf)/lw;                                              % Linear Mass of the cell

end