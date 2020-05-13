function [bw,T,status,xtraParams]=optThreshold(img,varargin)
% [bw,T,status]=optThreshold(img,varargin)
% perform thresholding of an image based on using defined algorithm
%
% Algorithm choice ('method'): 
%     {'minerr'},'otsu','kmeans','gm'
% 
% Param (not all are aplicable to all methods)
%
% 'lvl' => 2^12
% 'transform => {'none'} / 'log'
% 'solver' => 'fminsearch',{'fminunc'} (only relevent for some optimization
%                                       based methods)


%% define parames
arg.objsize=10;
arg.msk=true(size(img));
arg.method='minerr';
arg.transform='none';
arg.lvl=2^16; % some methods require discritization of the image, this param controls to how many levels
arg.solver='fminsearch';
arg.fallback='fminbnd';
arg.classnum=2;
arg.classnuminbackground=1;
arg.priorpercentforground=0.5;
arg.minnumofobjects=0;
arg.verbose=0;
arg.posteriorthresh=0.5;
arg.subsample=numel(img);
v=ver('stats');
arg.gmstart=[];
arg=parseVarargin(varargin,arg); 

%% init (make sure img is double etc.)
img=double(img);
bw=[];
%% first tranform the data as requested
switch arg.transform
    case 'none'
        gry=img(logical(arg.msk));
        gry=gry(:);
    case 'log'
        gry=log(img(img>0 & arg.msk));
    otherwise
        error('unsupported data tranform')
end

%% subsample (to make things faster) if required
if arg.subsample<numel(gry)
    gry=randsample(gry,arg.subsample);
end
    

%% deal with initial start for cases where method  gm
if strcmp(arg.method,'gm')
%     arg.gmstart = 'randSample';
    if str2double(v.Version) < 7.3 
        if ~isstruct(arg.gmstart) % make sure its a struct - must pass one into gmdistribution.fit
            arg.gmstart=struct(...
                'mu',prctile(gry(:),[25 75])',...
                'Sigma',cat(3,std(gry(gry<median(gry(:)))),std(gry(gry>median(gry(:))))),...
                'PComponents',[0.5 0.5]);
        end
    else 
        if ~isstruct(arg.gmstart) && numel(arg.gmstart) ~= numel(img) && arg.classnum==2
            if max(gry(:)>1) || min(gry(:))<0
                grytmp=mat2gray(gry);
            else
                grytmp=gry;
            end
            arg.gmstart = im2bw(grytmp,graythresh(grytmp));
            arg.gmstart=double(arg.gmstart(:))+1;
        end
    end
end

status=1;
arg.verbose && fprintf('started calculating threshold at %s\n',datestr(now)); %#ok<*VUNUS>
t0=now;

%% find threshold

switch arg.method
    case 'max'
        % include some optimization procedue that finds the cutof for
        % imextendedmax such that it maximizes the forground pixels num
    case 'robust'
        q=prctile(gry,[5 95]);
        indx=gry>q(1) & gry<q(2);
        T=mean(gry(indx))+2*std(gry(indx));
    case 'fminbnd'
        f=@(t) nnz(im2bw(img,prctile(img(arg.msk),t)) & ~bwareaopen(im2bw(img,prctile(img(arg.msk),t)),4));
        T=prctile(img(:),fminbnd(f,90,100));
%     case 'texture'
%         s=stdfilt(gry);
%         s2=backgroundSubtraction(s);
%         bw=optThreshold(s2);
%         bw=bw | im2bw(y,median(y(bw))); % combine std and intensity
%         bw=imclose(bw,clsSE);
%         bw=imfill(bw,'holes');
%         bw=bwareaopen(bw,mincellSize);
    case 'medmorph'
        % special case: use morphological operations as well
        T=prctile(gry,(1-arg.priorpercentforground)*100);
        bw = img>T;
        bw = bwareaopen(bw,arg.objsize,4);
        bw = imerode(bw,ones(2));
        bw = imclose(bw,ones(3));
        bw = bwareaopen(bw,arg.objsize,4);
        T=NaN;
        status=1;
    case 'edge'
        e=edge(img,'canny');
        bw = imdilate(e,ones(3));
        bw = imfill(bw,'holes');
        bw=imerode(bw,ones(3));
        T=NaN;
        status=1;
    case 'otsu'
        mn=min(gry);
        mx=max(gry);
        scl=mat2gray(gry);
        T=graythresh(scl)*(mx-mn)+mn; 
    case 'kmeans'
        [id,mu]=kmeans(gry(:),arg.classnum);
        [~,ordr]=sort(mu);
        bckid=ordr(arg.classnuminbackground);
        fgdid=ordr(arg.classnuminbackground+1);
        T=mean([min(gry(id==fgdid)); max(gry(id==bckid))]);
    case 'gm'
        try
            if ~isempty(arg.gmstart)
                gm=gmdistribution.fit(gry,arg.classnum,'start',arg.gmstart);
            else
                gm=gmdistribution.fit(gry,arg.classnum);
            end
            if nargout>3
                xtraParams=struct(gm);
            end
            lvl=linspace(min(gry(:)),max(gry(:)),arg.lvl);
            p=gm.posterior(lvl(:));
            [~,backgroundClass] = min(gm.mu);
            T=max(lvl(p(:,backgroundClass)>arg.posteriorthresh));
%             p=gm.posterior(lvl(:));
%             T=lvl(find(p(:,2)<p(:,1),1,'last'));
        catch ME 
           warning('gm failed - falling back to robust threshold, error was: %s',ME.message)
           q=prctile(gry,[5 95]);
           indx=img>q(1) & img<q(2);
           T=mean(gry(indx))+2*std(gry(indx));
        end
    case 'minerr'
        % create the histogram for the image
        bins=linspace(min(gry(:)),max(gry(:)),arg.lvl)';
        H=histc(gry(:),bins);
        H=H./sum(H);

        % initial guess (same as otsu)
        mn=min(gry);
        mx=max(gry);
        scl=mat2gray(gry);
        T0=graythresh(scl)*(mx-mn)+mn;
        
        switch arg.solver
            case 'fminunc'
                warning('not sure if fminunc works, need to specify gradient...'); %#ok<*WNTAG>
                if isempty(which('fminunc'))
                    error('solver not avaliable (do you have opt toolbox?) try running with fminsearch as solver instead');
                end
                opt=optimset('display','on','TolFun',1E-20);
                T = fminunc(@J,T0,opt);
            case 'fminsearch'
                [T,~,flag] = fminsearch(@J,T0);
                status=flag;
            case 'patternsearch'
                T = patternsearch(@J,T0);
            otherwise 
                error('unsupported solver');
        end
                
    otherwise
        error('unsupported thresholding method')
end

%% inverse transform threshold
switch arg.transform
    case 'none'
    case 'log'
        T=exp(T);
end
if isempty(bw)
    bw=img>T;
end
arg.verbose && fprintf('finish calculating threshold at %s total of %s\n',datestr(now),datestr(now-t0,13));
%% check to see if result are resonable and update status

% check to see if need to count objects
if arg.minnumofobjects>0
     [~,n]=bwlabel(bw);
end
if status~=true || (arg.minnumofobjects>0 && arg.minnumofobjects>n)
    warning('problem with thresholding - falling back to %s method',arg.fallback)
    [bw,T,status]=optThreshold(img,'method',arg.fallback,'msk',arg.msk,'transform',arg.transform);
end
        
            
% define the function to minimize 
function [y,grad]=J(T) %#ok<STOUT>
    ixb=bins<T;
    ixf=~ixb;
    Pb=sum(H(ixb));
    Pf=1-Pb;
    mub=1./Pb.*sum(H(ixb).*bins(ixb));
    muf=1/Pb*sum(H(ixf).*bins(ixf));
    sigb=sqrt(1/Pb*sum(H(ixb).*(bins(ixb)-mub)).^2);
    sigf=sqrt(1/Pf*sum(H(ixf).*(bins(ixf)-muf)).^2);
    y=1+2*(Pb*log(sigb)+Pf*log(sigf))-2*(Pb*log(Pb)+Pf*log(Pf));

end %of nested function

end % of main funciton