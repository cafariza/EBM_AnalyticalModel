%% Post-Processing / RESULT COMPARISON
% --------------------------------------------------------------------------

for i=1:length(T{1,1}(:,1))
    er_ex(i)=100*(abs(TCastem{1,1}(i,1))-abs(T{1,1}(i,1)))/TCastem{1,1}(i,1);
end
% Wall stiffness relation -Shear forces ??
%Te : Total shear force in nw walls /equivalent to Tr
TotSF_an=0;
for i=1:nw
TotSF_an=T{1}(1,i)+TotSF_an;
end 
TotSF_num=0;
for i=1:nw
TotSF_num=TCastem{1}(1,i)+TotSF_num;
end 