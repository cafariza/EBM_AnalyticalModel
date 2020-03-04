function [K,cmdout2] = ExCastemKM(K1,aw1,aw2,af1,af2,af,N,nw)
%This function helps to modify the text files that Castem uses to compute
%the shear stiffness of the story and the dynamic properties of a structure

t1=num2str(aw1);t2=num2str(aw2);t3=num2str(af1);t4=num2str(af2);t5=num2str(N);t6=num2str(nw);
text1='aw1=';text2='aw2=';text3='af1=';text4='af2=';tsep='; ';text5='N='; text6='nw=';
text_K=strcat(text1,t1,tsep,text2,t2,tsep,text3,t3,tsep,text4,t4,tsep,text6,t6,tsep);
text=strcat(text1,t1,tsep,text2,t2,tsep,text3,t3,tsep,text4,t4,tsep,text5,t5,tsep,text6,t6,tsep);

fileName2=sprintf('F2D_K_cbD.dgibi');
 S = fileread(fileName2);
S = [text_K, char(10), S];
FID = fopen(fileName2, 'w');
if FID == -1, error('Cannot open file %s', fileName2); end
fwrite(FID, S, 'char');
fclose(FID);
dir2=pwd;
pun1='CASTEM17 ';
pun='\';
filedir2 = strcat(pun1,{' '},dir2,pun,fileName2);
%command='CASTEM17 D:\Carolina\Thesis\Modeling\MatLab\F2D_K_cbD.dgibi';
command=filedir2{1};
[status,cmdout] = system(command);

fileName=sprintf('SHEARFORCE_b.inp');
startrow=1+nw*2+1+1+1;
endrow=startrow+2*nw-1;
ShearForcew=importK(fileName,startrow,endrow);
TSFi=zeros(nw*2,1);
TSF=0;
for i=1:nw*2
 SFi=textscan(ShearForcew{i}, '%u %f');
 TSFi(i)=SFi{2}*1E-6;
end
 TSFitemp=sort(TSFi);
for i=1:nw
  TSF=TSFitemp(i)+TSF;
end
K=TSF*-1;

fileName3=sprintf('1F_MOD_SY.dgibi');
S2 = fileread(fileName3);
S2 = [text, char(10), S2];
FID = fopen(fileName3, 'w');
if FID == -1, error('Cannot open file %s', fileName3); end
fwrite(FID, S2, 'char');
fclose(FID);
pun1='CASTEM17 ';
dir2=pwd;
pun='\';
filedir3 = char(strcat(pun1,{' '},dir2,pun,fileName3));

command=filedir3;
%command='CASTEM17 D:\Carolina\Thesis\Modeling\MatLab\1F_MOD_SY.dgibi';
[status,cmdout2] = system(command);

end

