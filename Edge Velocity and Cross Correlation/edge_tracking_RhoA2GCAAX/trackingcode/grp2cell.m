function [cl,unqid,avg]=grp2cell(M,id)
% function [cl,unqid,avg]=grp2cell(M,id)
% transform a grouped variable into cell array
% id must be numeric M a vector or a matrix
%
% as a "bonus" if nargout==3 it calcualtes the mean of ALL the elements in
% each cell (regardless is M is a ector aor a matrix!!! basically runs:
% f=@(x) mean(x(:))

if size(id,2)~=1, error('id must be a colum vector)'); end
if size(M,1)~=size(id,1), error('id must have same number of rows as M'); end

% remove id==0 
M=M(id>0,:);
id=id(id>0);

[id,ordr]=sort(id);
M=M(ordr,:); 
df=diff(id);
unqid=unique(id);
n=length(unqid);
df=find(df);
df=[0; df; length(id)];
cl=cell(n,1);
if size(M,2)==1 && exist('grp2cell_indxloop','file')==3
    if nargout==3
        [cl,avg]=grp2cell_indxloop(double(M),df);
    else
        cl=grp2cell_indxloop(double(M),df);
    end
else
    for i=1:n
        cl{i}=M((df(i)+1):df(i+1),:);
    end
    if nargout==3
        avg=cellfun(@(x) mean(x(:)),cl);
    end
end
    