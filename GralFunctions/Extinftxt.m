function [ MatrixData ] = Extinftxt( fileName, N,varargin)
%It extracts specific information from a text file 
if nargin<=2
nw=2;
end

nw=varargin{end};
startrow=1+N*nw+nw+1+1+1;
endrow=startrow+N*nw+nw-1;
data=importK(fileName,startrow,endrow);
dataco=importK(fileName,2,N*nw+nw+1);

MatrixData=zeros(N*nw+nw,3);
for i=1:N*nw+nw
 Coord=textscan(dataco{i}, '%f %f %f %f');
 Vdata=textscan(data{i}, '%f %f');
 MatrixData(i,:)=[(Coord{2}) (Coord{3}) abs(Vdata{2})];
end
MatrixData=sortrows(MatrixData,1);

end

