function [I,bck,bcksml] = backgroundSubtraction(I,varargin)
% background substruction of a segmented image. 
% 

sz = size(I);

if ~isa(I,'double') && ~isa(I,'single')
    warning('Data was not in double - tansforming it to be double using mat2gray!!!');
    I=mat2gray(I);
end

% default arguments
arg.smoothmethod='gauss'; % 'none' (note that smoothing also take cares of NaN if exist)  - used to be tpaps
arg.msk=ones(sz); % other can be any logical mask of what could be background
arg.samplingmethod='block'; % 'grid' 
arg.samplingfcn=@(x) min(x.data(:));%@(x) prctile(x(:),5); %@(x) min(x(:)) / @(x) median(x(:))
arg.samplingdensity=15;
arg.fgauss=fspecial('gauss',11,5);
arg.interpmethod='bilinear';

arg=parseVarargin(varargin,arg);
% for i=1:2:length(varargin)
%     arg.(varargin{i})=varargin{i+1};
% end

Imsk=I; 
Imsk(~arg.msk)=NaN;

%% sample bckgroud pixels using either grid or blocks (~mask = NaN)
switch arg.samplingmethod
    case 'grid'
        % set up sampling grid. 
        c = floor(linspace(1,sz(2),arg.samplingdensity));
        r = floor(linspace(1,sz(1),arg.samplingdensity));
        [c,r]=meshgrid(c,r);

        % sample from block within the image
        bcksml=zeros(size(r));
        bcksml(:)=Imsk(sub2ind(sz,r(:),c(:)));
    case 'block'
        blksz=bestblk(sz,max(sz)/arg.samplingdensity);
        bcksml=blockproc(Imsk,blksz,arg.samplingfcn);
        r=(blksz(1)/2):blksz(1):sz(1) ;
        c=(blksz(2)/2):blksz(2):sz(2);
        [c,r]=meshgrid(c,r);
    otherwise
        error('unsupported samplind method')
end

%% Smooth (also removes some of the NaNs along the way)
ix=find(~isnan(bcksml));
switch arg.smoothmethod
    case 'none'
        msksml=imresize(arg.msk,size(bcksml));
        bcksml(~msksml)=mean(I(arg.msk>0));
    case 'griddata'
        % fit a
        [x,y]=meshgrid(1:sz(2),1:sz(1));
        bcksml=griddata(c(ix),r(ix),bcksml(ix),x,y);
        % get rid of the NaNs at the edge of the griddata interpolation
        [rstart,cstart]=find(~isnan(bcksml),1);
        [rend,cend]=find(~isnan(bcksml),1,'last');
        bcksml=bcksml(rstart:rend,cstart:cend);
        bcksml=imresize(bcksml,sz);
    case 'triscatter'
        F = TriScatteredInterp(c(ix),r(ix),bcksml(ix));
        [x,y]=meshgrid(1:sz(2),1:sz(1));
        bcksml = F(x,y);
    case 'spline'
        % thin plate spline fitting to points
        st = tpaps([c(ix) r(ix)]',bcksml(ix)');
        avals = fnval(st,[c(:) r(:)]'); % evaluate tps function
        % transform tps values from vector to matrix form
        bcksml = reshape(avals,size(r));
    case 'gauss'
        bcksml = imfilter(bcksml,arg.fgauss,'symmetric');
    otherwise
        error('unsupported smoothing method')
end

%% Interpolate back to big size
bck=imresize(bcksml,sz,arg.interpmethod);
I=imsubtract(I,bck);
I(I<0)=0;


