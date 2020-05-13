function [ argstruct ] = parseVarargin( args,argstruct )
if ~isempty(args)
    for i=1:2:length(args)
    argstruct.(args{i})=args{i+1};
    end
end
