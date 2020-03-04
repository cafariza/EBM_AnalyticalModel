function [ NewMatrixData ] = ExtinftxtForces( fileName, N)
%It extracts specific information from a text file of the form of CASTEM
%FORCES Example: IFN1MOD1.inp

fid = fopen(fileName);
res={};
while ~feof(fid)
  thisline = fgetl(fid);
  if ~ischar(thisline); break; end
  res{end+1,1} = thisline;
end
fclose(fid);
Numblines=numel(res);
Nwalls=(Numblines-8)/(2*(N+1));
startrow=1+N*Nwalls+Nwalls+1+7;
endrow=startrow+N*Nwalls+Nwalls-1;
data=importK(fileName,startrow,endrow);
dataco=importK(fileName,2,N*Nwalls+Nwalls+1); %OK
MatrixData=zeros(N*Nwalls+Nwalls,5);  %OK
for i=1:(N*Nwalls+Nwalls)
 Coord=textscan(dataco{i}, '%f %f %f %f');
 Vdata=textscan(data{i}, '%f %f %f %f %f %f %f');
 MatrixData(i,:)=[(Coord{2}) (Coord{3}) (Vdata{2}) (Vdata{3}) (Vdata{7})];
end
MatrixData=sortrows(MatrixData,1);
%Splitting
NewMatrixData={};

startrow=1;

for i=1:Nwalls
endrow=startrow+N;
NewMatrixData{i}=MatrixData(startrow:endrow,:);
startrow=endrow+1;
for j=1:N+1
 if j==1
      NewMatrixData{i}(j,:)=NewMatrixData{i}(j,:);
 else
      NewMatrixData{i}(j,:)=[NewMatrixData{i}(j,1) NewMatrixData{i}(j,2) (NewMatrixData{i}(j,3)*2-NewMatrixData{i}(j-1,3)) (NewMatrixData{i}(j,4)*2-NewMatrixData{i}(j-1,4)) (NewMatrixData{i}(j,5)*2-NewMatrixData{i}(j-1,5))];
 end
end
end
end

