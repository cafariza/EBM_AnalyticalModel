function [K,cmdout2] = ExCasteminf(fileName2,Iy,Iz,Sect1,poisson_eq,N_elem)
%This function helps to modify the text files that Castem uses to compute
%the shear stiffness of the story and the dynamic properties of a structure
t1=num2str(N_elem);t2=num2str(Sect1/10^5);t3=num2str(Iz/10^5);t4=num2str(Iy/10^5);t5=num2str(poisson_eq);
text1='NbN=';text2='Sect1=';text3='Iz=';text4='Iy=';tsep='; ';text5='Nustru='; E5='E5'; 
text_K=strcat(text1,t1,tsep,text2,t2,E5,tsep,text3,t3,E5,tsep,text4,t4,E5,tsep,text5,t5,tsep);
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
command=filedir2{1};
[status,cmdout] = system(command)
end