function []=savecase(aw1,af1,N,nw)
t1=num2str(aw1);
t3=num2str(af1);
t5=num2str(N);
t6=num2str(nw);
Rootdir = pwd;


nextdir= (['\N' t5 'nw' t6 t1 t3]);
destdirectory= ([Rootdir '\MainvarFEMCasesNovember2019' nextdir ]);
mkdir(destdirectory);   %create the directory
for k=1:13
baseFileName = sprintf('figure_%d.jpg',k);
baseFileNamefig = sprintf('figure_%d.fig',k);
% Specify some particular, specific folder:
fullFileName = fullfile(destdirectory, baseFileName); 
fullFileNamefig = fullfile(destdirectory, baseFileNamefig);
saveas(figure(k),fullFileName); % Using export_fig instead of saveas.
saveas(figure(k),fullFileNamefig); % Using export_fig instead of saveas.
end
end 