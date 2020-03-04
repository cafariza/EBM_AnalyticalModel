Function 
t1=num2str(aw1);t2=num2str(aw2);t3=num2str(af1);t4=num2str(af2);t5=num2str(N);t6=num2str(nw);
text1='aw1=';text2='aw2=';text3='af1=';text4='af2=';tsep='; ';text5='N='; text6='nw=';
text_K=strcat(text1,t1,tsep,text2,t2,tsep,text3,t3,tsep,text4,t4,tsep);
text=strcat(text1,t1,tsep,text2,t2,tsep,text3,t3,tsep,text5,t5,tsep,text6,t6,tsep);

if length(af)<3
    K
else
fileName2=sprintf('FRAME2D_K_cbD.dgibi');
S = fileread(fileName2);
S = [text_K, char(10), S];
FID = fopen(fileName2, 'w');
if FID == -1, error('Cannot open file %s', fileName2); end
fwrite(FID, S, 'char');
fclose(FID);

command='CASTEM17 D:\Carolina\Thesis\Modeling\Codes MatLab\GeneralFunctions\FRAME2D_K_cbD.dgibi';
[status,cmdout] = system(command);

fileName=sprintf('SHEARFORCE_b.inp');
startrow=1+length(aw)*2+1+1+1;
endrow=startrow+length(aw)-1;
ShearForcew=importK(fileName,startrow,endrow);
TSFi=0;
for i=1:length(aw)
 SFi=textscan(ShearForcew{i}, '%u %f');
 TSFi=SFi{2}*1E-6+TSFi;
end
K=TSFi*-1
end

fileName3=sprintf('1F_MODAL.dgibi');
S2 = fileread(fileName3);
S2 = [text, char(10), S2];
FID = fopen(fileName3, 'w');
if FID == -1, error('Cannot open file %s', fileName3); end
fwrite(FID, S2, 'char');
fclose(FID);
command='CASTEM17 D:\Carolina\Thesis\Modeling\Codes MatLab\1F_MODAL.dgibi';
[status,cmdout] = system(command)