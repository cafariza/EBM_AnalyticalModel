function []=savecase(aw1,af1,N,nw)
%Parameters for saved images
ImageDPI=500;
ImageSizeX=7;
ImageSizeY=5;
t1=num2str(aw1);
t3=num2str(af1);
t5=num2str(N);
t6=num2str(nw);
Rootdir = pwd;


nextdir= (['\N' t5 'nw' t6 t1 t3]);
destdirectory= ([Rootdir '\MainvarFEMCasesNovember2019' nextdir ]);
mkdir(destdirectory);   %create the directory
figs = findobj(0, 'type', 'figure'); 
for k=1:length(figs)
FileLabel=sprintf('figure_%d',k);
baseFileNamefig = sprintf('figure_%d.fig',k);
% Specify some particular, specific folder:
fullFileNamefig = fullfile(destdirectory, baseFileNamefig);
saveas(figure(k),fullFileNamefig); % Using export_fig instead of saveas.
%Set the image size and save to FileLabel.png where FileLabel is set at line 9. 
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 ImageSizeX ImageSizeY])
print('-dpng', strcat(FileLabel, '.png') , strcat('-r',num2str(ImageDPI)))
end
%====================

end 