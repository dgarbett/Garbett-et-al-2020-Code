function opt=setopts(opt,varargin)

if mod(length(varargin),2)~=0
    error('options should always be in pairs')
end
for i=1:2:length(varargin)
    o=varargin{i}; % o is what optio to change
    v=varargin{i+1}; % v is what is its value
    
    % o must be a char
    if ~ischar(o)
        error('Option must be a char!');
    end
    
    % check if option is an allowed one
    if ~ismember(o,fieldnames(opt))
        error('%s is not a valid option!',o)
    end
      
    opt.(o)=v;
end
