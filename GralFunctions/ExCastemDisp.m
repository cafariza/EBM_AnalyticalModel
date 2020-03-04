function [UY_numuA,UY_numuU] = ExCastemDisp(nw)
%This function helps to modify the text files that Castem uses to compute
% 

fileName=sprintf('DISPY_u.inp');
startrow=1+nw*2+1+1+1;
endrow=startrow+nw-1;
DISPY_u=importK(fileName,startrow,endrow);
Coord=importK(fileName,2,2+nw-1);
infoth=zeros(nw,3);
for i=1:nw
 Vdisp=textscan(DISPY_u{i}, '%f %f');
 Coordx=textscan(Coord{i}, '%f %f %f %f');
 infoth(i,1:end)=[Vdisp{1} Coordx{2} Vdisp{2}];
end
infotho=sortrows(infoth,2);
UY_numuU=infotho(:,3);

fileName=sprintf('DISPY_a.inp');
DISPY_a=importK(fileName,startrow,endrow);
Coord=importK(fileName,2,2+nw-1);
infoth=zeros(nw,3);
for i=1:nw
 Vdisp=textscan(DISPY_a{i}, '%f %f');
 Coordx=textscan(Coord{i}, '%f %f %f %f');
 infoth(i,1:end)=[Vdisp{1} Coordx{2} Vdisp{2}];
end
infotho=sortrows(infoth,2);
UY_numuA=infotho(:,3);

end

