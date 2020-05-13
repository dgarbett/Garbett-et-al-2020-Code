function filenames=getFilenames(path)
% returns a list of the filenames in a specified directory

d0=dir(path);
filenames={d0([d0.isdir]==0).name};
filenames=filenames(:);
filenames=setdiff(filenames,{'Thumbs.db','thumbs.db'});
