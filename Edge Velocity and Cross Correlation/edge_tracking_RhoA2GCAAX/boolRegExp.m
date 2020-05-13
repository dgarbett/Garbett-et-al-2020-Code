function bool=boolRegExp(cell_arr,pattern)
%Returns a boolean vector or matrix indicating whether or not the elements
%of cell_arr contain pattern

bool=cellfun(@(x) length(x),regexp(cell_arr,pattern))>0;