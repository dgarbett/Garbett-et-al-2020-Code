function M=annsearchsingle(Xr, Xq, maxdisp)
% wraps annsearch but returns a little different output
% M has the query indx and reference index and distance

% This is an adaptation of annquerysingle intended to work with annsearch
% instead of annquery

%[ix,d]=annsearch(Xr,Xq,1,'search_sch','fr','radius', maxdisp);
[ix,d]=annsearch(Xr,Xq,1);
N=1:size(Xq,2);
ix=single(ix(d<=maxdisp))';
N=single(N(d<=maxdisp))';
d=single(d(d<=maxdisp))';
M=sortrows([N ix d]);
if isempty(M)
    return
end
[bla,ind]=unique(M(:,2),'first');
M=M(ind,:);
