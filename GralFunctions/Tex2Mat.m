function [Mat]= Tex2Mat(f_path,file_t1,file_t3,n,m)
% Matrices generator from text file
for i=1:length(n)
    file_t2=num2str(n(i));
    fname=strcat(f_path,file_t2);
    path ([pwd fname],path);
    for j=1:length(m)
        file_t4=num2str(m(j));
        fid=fopen(strcat(file_t1,file_t2,file_t3,file_t4));
        tLines = fgets(fid);
        numCols = (numel(strfind(tLines,';')))+(numel(strfind(tLines,',')))+12;
        C=textscan( fid, repmat('%f',[1,numCols*2+1]),1 ...
                ,   'CollectOutput' ,   true    ...
                ,   'Delimiter'     , ',;'     );
        
        fclose(fid);  
        fid=fopen(strcat(file_t1,file_t2,file_t3,file_t4));
                C=textscan( fid, repmat('%f',[1,numCols*2+1]) ...
                ,   'CollectOutput' ,   true    ...
                ,   'Delimiter'     , ',;'     );
        
        fclose(fid);  
        m_data=C{1,1};
        s1=size(m_data,1);
        for k= 1: s1
          msubs=m_data(k,1:end); % Save existing data in ith row of m_data
            msubs=msubs(isnan(m_data(k,1:end))==0); %Substitute matrix/ taking only non-NaN values
            m_data(k,1:end)=0; %Erase all existing values in ith row of m_data
            m_data(k,1:size(msubs,2))=msubs; %Substitute values without NaN
        end
        Mat{i,j}=m_data; %Cell array with all the data info
    end    
end
save Mat
end