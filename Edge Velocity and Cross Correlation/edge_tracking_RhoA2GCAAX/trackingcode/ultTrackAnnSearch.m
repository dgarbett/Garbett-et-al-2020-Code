function xyTraj = ultTrackAnnSearch(pnt,varargin)
% Create trajectories from xy points in multiple frames. 
% pnt is a cell array such that each cell is the points in one frame
% in each frame points are XYZI etc. 
% 
% optional arguments: 
% 
% T   => a vector of the time for each frame.
%        if not provided, default is {1:length(pnt)}
% dim => which colum in the pnt matrices to use for tracking
% mem => how many (if at all) frames to  "skip" and connect traj {2}
% fmemmaxdist => a function handle that gets mem and disp and retun the
%                max distance to consider for that mem 
%                {@(mem,mxdist) sqrt(mem)*mxdist}
% maxdisp => max displacment: a number {5}
% verbose = > verbode: {true} / false
% jitter => jitter correction method (could be more than one in a cell
% array:
%      {'none'}       no correction what so ever
%       'meanDrift'   fits a spline to mean xy positions
%       'medianDrift' fits a splint to median xy
%       'meanFull'    make sure the median of all cells is the same
%       'medianFull'  make sure the median of all cells is the same
%                                 
% minlength => portion of total number of frames need to be a legit raj {0} (goes from 0-1)
% pairrule => what rule to use for pairings: 
%       {'fwdnn'}     is the forward nereaset neighbour  
%        'bcknn'      is a two background NN
%        'fwdbckmtch' is a two directional NN, e.g. cell in frame i
%                     is the nn to cell in frame j and also vice
%                     versa. Only pairs is bck and fwd agree.
%        'maxflow'    uses max-flow algorithm for the mattching,
%                     its nice but 1. very slow (O(N^2*E)) and requires
%                      bioinformatics toolbox
%        'seednn'     track only in reference to a seed set of points.
%        'hungarian'  Uses the hungarian algorithm to match points between
%                     frames (uses the hungarianAssignment function). It is
%                     a non-greedy algorithm. It works only on all the
%                     pairs of points that are maxdisp or smaler from one
%                     another. 


% Written by Roy Wollman, Oct 2008. 

%% deal with input arguemnts
t0=now;

% define defaults for all input arguments
arg.dim=1:2; % dim => which colum in the pnt matrices to use for tracking
arg.mem=0; % how many (if at all) frames to allow a "skip" and connect trajectory
arg.maxdisp=5; % max displacment: a number {5}
arg.verbose=true; % verbode: {true} / false
arg.jitter='none'; % jitter correction: true / {false}
arg.minlength=0; % a number (0-1) of the portion of the number of frames a traj eed to be to be considered legit
arg.pairrule='fwdnn'; % chose of pairing algorithm
arg.maxframe=Inf;
arg.T=1:length(pnt);
arg.seed=[];
arg.fmemmaxdist=@(mem,mxdist) sqrt(mem)*mxdist;

fld=fieldnames(arg);
for i=1:2:length(varargin)
    if ~any(ismember(fld,varargin(i)))
        error('non supported input argument'); 
    end
    arg.(varargin{i})=varargin{i+1};
end

% if pnt is char interpet it as a path and read point from file
if ischar(pnt)
    arg.verbose && fprintf('Started reading data from disk %s\n',datestr(now-t0,13));
    [pnt,arg.T]=readCellPointData(pnt); 
    arg.verbose && fprintf('Finished reading data from disk %s\n',datestr(now-t0,13));
end

% truncate pnt to max frame
pnt=pnt(1:min(length(pnt),arg.maxframe));
pnt=pnt(:).';            

%% calculate size vecotrs and preallocate memory
n=cellfun(@(x) size(x,1),pnt);
%n=cellfun(@length,pnt);
ncmsm=cumsum([0 n]);
ncmsm=ncmsm(:);

frameNum=length(pnt);
match=cell(frameNum-1,1);

%% jitter correction
% multiple possible methods, can run more than one, 
if ~iscell(arg.jitter), arg.jitter={arg.jitter}; end
for i=1:length(arg.jitter)
    switch arg.jitter{i}
        case 'none'
            jitterCorrection=[];
        case {'meanDrift','medianDrift'}
            switch arg.jitter{i}
                case 'meanDrift'
                    mn=cellfun(@(x) mean(x(:,arg.dim)), pnt,'uniformoutput',false);
                case 'medianDrift'
                    mn=cellfun(@(x) mean(x(:,arg.dim)), pnt,'uniformoutput',false);
            end
            mn=cat(1,mn{:});
            ft=zeros(size(mn,1),length(arg.dim));
            for j=1:length(arg.dim)
                spln = spap2(1,4,1:size(mn,1),mn(:,arg.dim(j)));
                ft(:,j)=fnval(spln,1:size(mn,1));
            end
            jitterCorrection = [0 0; diff(ft)];
            jitterCorrection = cumsum(jitterCorrection);
        case {'meanFull','medianDrift'}
            switch arg.jitter{i}
                case 'meanFull'
                    mn=cellfun(@(x) mean(x(:,arg.dim)), pnt,'uniformoutput',false);
                case 'medianFull'
                    mn=cellfun(@(x) mean(x(:,arg.dim)), pnt,'uniformoutput',false);
            end
            jitterCorrection = [0 0; diff(mn)];
            jitterCorrection = cumsum(jitterCorrection);
        otherwise
            error('unsopported jitter correction method')
    end
end

for i=1:size(jitterCorrection,1)
    pnt{i}(:,arg.dim)=pnt{i}(:,arg.dim)-repmat(jitterCorrection(i,:),n(i),1);
end

%% create pairing between frames using specified algorithm 

arg.verbose && fprintf('Started calculation of pairs at %s\n',datestr(now-t0,13));

switch arg.pairrule
    case 'seednn'
        % first find in the first frame the neaerest neighbour poitns to
        % the seeds
        M=annsearchsingle(double(arg.seed'),double(pnt{1}(:,arg.dim)'),arg.maxdisp);
        match{1}(:,1)=double(M(:,1));
        for i=2:frameNum
            M=annsearchsingle(double(arg.seed'),double(pnt{i}(:,arg.dim)'),arg.maxdisp);
            mtch=1:length(arg.seed);
            mtch(M(:,2))=M(:,1)+(i-1)*length(arg.seed);
            match{i-1}(:,2)=mtch;
            match{i}(:,1)=mtch;
        end
    case 'fwdnn'
        for i=1:frameNum-1
            % find based on ann
            M=annsearchsingle(double(pnt{i+1}(:,arg.dim)'),double(pnt{i}(:,arg.dim)'),arg.maxdisp);
            match{i}=double(M(:,1:2))+repmat(ncmsm(i:i+1)',size(M,1),1);
        end
    case 'bcknn'
        for i=1:frameNum-1
            % find based on ann
            M=annsearchsingle(pnt{i}(:,arg.dim)',pnt{i+1}(:,arg.dim)',arg.maxdisp);
            match{i}=double(M(:,[2 1]))+repmat(ncmsm(i:i+1)',size(M,1),1);
        end
    case 'fwdbckmtch'
        for i=1:frameNum-1
            % find based on ann
            M1=annsearchsingle(pnt{i+1}(:,arg.dim)',pnt{i}(:,arg.dim)',arg.maxdisp);
            M2=annsearchsingle(pnt{i}(:,arg.dim)',pnt{i+1}(:,arg.dim)',arg.maxdisp);
            M=intersect(M1(:,1:2),M2(:,[2 1]),'rows');
            match{i}=double(M(:,1:2))+repmat(ncmsm(i:i+1)',size(M,1),1);
        end
    case  'hungarian'
        for i=1:frameNum-1
            [N,D]=annsearch(pnt{i+1}(:,arg.dim)',pnt{i}(:,arg.dim)',size(pnt{i+1},1),'search_sch','fr','radius', arg.maxdisp);
            CostMat=inf(size(pnt{i},1),size(pnt{i+1},1));
            for j=1:size(N,2)
                CostMat(j,N(N(:,j)>0,j))=D(N(:,j)>0,j)';
            end
            CostMat2=CostMat.*CostMat;   % added by Sean on 4-25-09 to make the cost be the squared distance instead of the distance
            assignment = hungarianAssigement(CostMat2);  %changed to CostMat2 by Sean
            id=1:size(pnt{i},1);
            ix=find(assignment>0);
            id=id(ix);
            assignment=assignment(ix);
            d=CostMat(sub2ind(size(CostMat),id,assignment));
            ix=find(d < arg.maxdisp);
            id=id(ix);
            assignment=assignment(ix);
            M=[id(:) assignment(:)];
            match{i}=double(M(:,1:2))+repmat(ncmsm(i:i+1)',size(M,1),1);
        end
    case 'maxflow'
        error('max flow not yet supported');
        for i=1:frameNum-1 %#ok<UNRCH>
            tri1=delaunayn(pnt{i}(:,arg.dim));
            [ix1,d1]=dsearchn(pnt{i}(:,arg.dim),tri1,pnt{i+1}(:,arg.dim));
            tri2=delaunayn(pnt{i+1}(:,1:2));
            [ix2,d2]=dsearchn(pnt{i+1}(:,arg.dim),tri2,pnt{i}(:,arg.dim));

            % remove nn that are further than arg.maxdisp away
            N1=1:n(i);
            N2=n(i)+1:n(i)+n(i+1);
            ix1=ix1(d1<=arg.maxdisp);
            N2=N2(d1<=arg.maxdisp);
            ix2=ix2(d2<=arg.maxdisp);
            N1=N1(d2<=arg.maxdisp);

            % create the graph
            % 1:n(i) are cells in frame i
            % n(i)+1:n(i+1) are cells in frame i+1
            % n(i)+n(i+1)+1 is the source
            % n(i)+n(i+1)+2 is the target for the max flow
            gn=n(i)+n(i+1)+2; % graph size
            G=sparse(gn,gn);
            % nearst neighbours of cells in frame i+1 to cell in frame i
            G(sub2ind([gn gn],N1,n(i)+ix2'))=1;
            % nearst neighbours of cells in frame i to cell in frame i+1
            G(sub2ind([gn gn],ix1',N2))=1;
            % source to cell in i
            G(gn-1,1:n(i))=1;
            G(gn,n(i)+1:n(i)+n(i+1))=1;
        end
    otherwise
        error('Not a sported pairing algorithm supported');
end
arg.verbose && fprintf('Fisniehd calculation fo pairs at %s\n',datestr(now-t0,13));

%% connect pairs into trajectories
% I'm using the match point to create a pairing graph implemented as a vector.
% the index is left node and the value is the right node. Than I'm just
% traversing this graph assigning id to each trajectory (the mex function)

arg.verbose && fprintf('Linking pairs into trajectories at %s \n',datestr(now-t0,13));
arg.verbose && fprintf('creating graph %s \n',datestr(now-t0,13));

gr=uint32(zeros(1,ncmsm(end)));
possibleStart=false(1,ncmsm(end));
match=uint32(cat(1,match{:}));
gr(match(:,1))=match(:,2);
possibleStart(match(:,1))=true;

% clear memory as we go
clear match

% link pairs into trajectories by traversing the graph
arg.verbose && fprintf('walking graph to create trajectories %s \n',datestr(now-t0,13));

if exist('traverseGraph','file')==3
    % make sure that gr and possibleStart are row vectors
    gr=gr(:).';
    possibleStart=possibleStart(:).';
    traj=traverseGraph(gr,possibleStart);
else % do it without the mex file
    cnt=0;
    traj=zeros(ncmsm(end),1);
    while any(possibleStart)
        cnt=cnt+1;
        nxt=find(notvisited,1,'first');
        while nxt>0
            traj(nxt)=cnt;
            possibleStart(nxt)=false;
            nxt=gr(nxt);
        end
    end
end

% clear memory as we go
clear gr possibleStart

%% convert traj map into cell array
nds=1:ncmsm(end);
nds=nds(traj>0);
traj=traj(traj>0);
allTraj=grp2cell(nds(:),traj(:));
arg.verbose && fprintf('Finished linking at %s \n',datestr(now-t0,13));

%% connect broken trajectories 
% logic of this section: first create the xy position of all the ends and
% starts of trajectories. Than go over all finish and start to start+mem
% and see if they are worth connectin. use same nearest neighbour rule as
% initial connection

% create a cew representation of trajectories that has all data from pnt 
% plus traj number and frame number. 

% note that for frame number we're using the trjStart + seq afterwardss and
% it works since at this point there were no "memory" jumps. 
XY=cat(1,pnt{:});
%trjStart=cellfun(@(traj) find(traj(1)<ncmsm,1)-1,allTraj);
trjStart=cellfun(@(traj) find(traj(1)<=ncmsm,1)-1,allTraj);
xyTraj = cellfun(@(x,y,z) [XY(x,:) z*single(ones(size(x,1),1)) single((y:y+size(x,1)-1))'],...
                   allTraj,... the trajectory ids
                   num2cell(trjStart),... the start frame to use for frame numbers
                   num2cell((1:length(allTraj))'),... trajectory ID
                   'uniformoutput',false);

if arg.mem==0, return, end
arg.verbose && fprintf('Connecting broken trajectories at %s \n',datestr(now-t0,13));

% create the start/finish point of trajectories
strt=cellfun(@(x) x(1,[arg.dim end-1:end]),xyTraj,'uniformoutput',false);
strt=cat(1,strt{:});
fnsh=cellfun(@(x) x(end,[arg.dim end-1:end]),xyTraj,'uniformoutput',false);
fnsh=cat(1,fnsh{:});
[strtcl,strtfrm]=grp2cell(strt,strt(:,end));
[fnshcl,fnshfrm]=grp2cell(fnsh,fnsh(:,end));

% create the list of trajectories that could be connected 
% based on startfnish points arg.mem and arg.maxdisp
% what this loop is doing is creating a matrix connectTraj with all
% trajectories that are 'legal' to connectwith the folosing form: 
% [fnshTrajIdx strtTrajIdx distance dframe];
connectTraj=cell(frameNum-1,1);
for i=1:frameNum-1
    % find all the xy positions for points that and at fnsdhfrm and start
    % at 2-mem afterwords
    xyf=cat(1,fnshcl{fnshfrm==i});
    xys=cat(1,strtcl{ismember(strtfrm,i+1+(1:max(1,arg.mem-1)))});
    if isempty(xyf) || isempty(xys)
        continue
    end
    M=annsearchsingle(double(xys(:,arg.dim)'),double(xyf(:,arg.dim)'),arg.fmemmaxdist(arg.mem,arg.maxdisp));
    if ~isempty(M)
        connectTraj{i}=[xyf(M(:,1),3) xys(M(:,2),3) M(:,3) xys(M(:,2),4)-xyf(M(:,1),4)];
    end
end
connectTraj=cat(1,connectTraj{:});
arg.verbose && fprintf('found connections, recreating xyTraj at %s \n',datestr(now-t0,13));

% Choose which connection to perform (e.g. if two traj want to link to the
% same one, choose the one that has the smaller frame jump and closest by distance)

% sort first by dFrame than by distance
connectTraj=sortrows(connectTraj,[4 3]);

% I can only allow single connection per start
[bla,idx]=unique(connectTraj(:,2),'first');
connectTraj=connectTraj(idx,:);

if ~isempty(connectTraj)
    % use traverseGraph again to identify which traj should be togather
    gr=uint32(zeros(1,length(xyTraj)));
    possibleStart=true(1,length(xyTraj));
    connectTraj=uint32(connectTraj(:,1:2));
    gr(connectTraj(:,1))=connectTraj(:,2);
    gr=gr(:).';
    possibleStart=possibleStart(:).';
    traj2connect=traverseGraph(gr,possibleStart);
    nds=1:length(allTraj);
    traj2connect=grp2cell(nds(:),traj2connect(:));
    finalAllTraj=cell(length(traj2connect),1);
    for i=1:length(traj2connect)
        finalAllTraj{i}=cat(1,allTraj{traj2connect{i}});
    end
end

% recalculate the xyTraj to return it in a reasonable way. 
% trjStart=cellfun(@(traj) find(traj(1)<ncmsm,1)-1,allTraj);
% xyTraj = cellfun(@(x,y) [XY(x,:) (y:y+size(x,1)-1)'],...
%                            allTraj,... the trajectory ids
%                            num2cell(trjStart),... the start frame to use for frame numbers
%                            'uniformoutput',false);

arg.verbose && fprintf('calculatig time vectors: %s\n',datestr(now-t0,13));

% expand the arg.T to have an element for each poin in the nerwork 
arg.T=arg.T(1:length(n));
TT=arrayfun(@(n,i) i*ones(n,1),n(:),arg.T(:),'uniformoutput',false); 
TT=cat(1,TT{:});
% find for each traj what its time points are. 
T=cellfun(@(x) TT(x),finalAllTraj,'uniformoutput',false);

% TT=[]; 
% for i=1:length(n), 
%     TT=[TT; arg.T(i)*ones(n,1)]; 
% end

% c = cellfun(@colon,num2cell(start),num2cell(stop),'UniformOutput',0);
% x = cat(2,c{:});

xyTraj = cellfun(@(x,y) [XY(x,:) y(:)],finalAllTraj,T,'uniformoutput',false);


arg.verbose && fprintf('done %s \n',datestr(now-t0,13));
