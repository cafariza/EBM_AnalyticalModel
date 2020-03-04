function [Theta_numuA, Theta_numuU] = ExCastemTheta(aw1,aw2,af1,af2,nw)
%This function helps to modify the text files that Castem uses to compute
% THETA
t1=num2str(aw1);t2=num2str(aw2);t3=num2str(af1);t4=num2str(af2);t6=num2str(nw);
text1='aw1=';text2='aw2=';text3='af1=';text4='af2=';tsep='; '; text6='nw=';
text_theta=strcat(text1,t1,tsep,text2,t2,tsep,text3,t3,tsep,text4,t4,tsep,text6,t6,tsep);

fileName2=sprintf('F2D_U_cbD_dens.dgibi');
S = fileread(fileName2);
S = [text_theta, char(10), S];
FID = fopen(fileName2, 'w');
if FID == -1, error('Cannot open file %s', fileName2); end
fwrite(FID, S, 'char');
fclose(FID);
dir2=pwd;
pun1='CASTEM17 ';
pun='\';
filedir2 = strcat(pun1,{' '},dir2,pun,fileName2);
%command='CASTEM17 D:\Carolina\Thesis\Modeling\MatLab\F2D_U_cbD_dens.dgibi';
command=filedir2{1};
[status,cmdout] = system(command);

fileName=sprintf('Theta_U.inp');
startrow=1+nw*2+1+1+1;
endrow=startrow+nw-1;
Theta_U=importK(fileName,startrow,endrow);
Coord=importK(fileName,2,2+nw-1);
infoth=zeros(nw,3);
for i=1:nw
 Vtheta=textscan(Theta_U{i}, '%f %f');
 Coordx=textscan(Coord{i}, '%f %f %f %f');
 infoth(i,1:end)=[Vtheta{1} Coordx{2} Vtheta{2}];
end
infotho=sortrows(infoth,2);
Theta_numuU=infotho(:,3);


fileName=sprintf('Theta_a.inp');
startrow=1+nw*2+1+1+1;
endrow=startrow+nw-1;
Theta_a=importK(fileName,startrow,endrow);
Coord=importK(fileName,2,2+nw-1);
infoth=zeros(nw,3);
for i=1:nw
 Vtheta=textscan(Theta_a{i}, '%f %f');
 Coordx=textscan(Coord{i}, '%f %f %f %f');
 infoth(i,1:end)=[Vtheta{1} Coordx{2} Vtheta{2}];
end
infotho=sortrows(infoth,2);
Theta_numuA=infotho(:,3);

end

