%% Plate_single_cell_speed_4x_just leaders_Yilin-whist.m
%adapted from Arnold's code with help(lots) from Yilin to find edge of
%scratches and track only cells within set number of pixels from edge (100
%works well)
%% init parameters *CHANGE for each experiment*
clear; close all
clc;
% folder where data is located
root='E:\Paper Data\Raw\siRNA scratch assay\migration rates\06.16.16\20160616-siRNA\Results\';
rawdir= [root filesep];
filenames=getFilenames(rawdir); %gets tracedata filenames in root directory
datadir=[root,'Strip100_TP3-70_Analysis_' date];
mkdir(datadir);
%% acquisition params *CHANGE for each experiment*
st=3; % start frame
ed=70; % end frame
pix=1.62; %IXMicro 20x no bin=0.32, 20x 2x2= 0.65, 10x no bin=0.65, 10x 2x2=1.31, 4x no bin=1.62
tstep=1; % time step (1= movement at every tp, 2= movement every other tp)
strip=100; % width of bands of cells from sheet boundaries to be analyzed, ie... how deep into sheet from edge do we analyze (in pixels)?

%% goes through tracedata files in root directory (files are previously generated from Timelapse.m)

for n=1:length(filenames)
        tic
        tempfile=filenames{n,1}; %load filename
        table(n).filename=tempfile;
        load([rawdir tempfile]); %loads data from file
        lo=find(~isnan(tracedata(:,st,1))&~isnan(tracedata(:,ed,1))); % finds tracks going from st to ed
        
        Y=tracedata(lo,st,2);% Y-values of cells at ed.  %X=tracedata(lo,st,1);
        sep=kmeans(Y,2,'start',[1;max(Y)]); % seperates all cells in two groups according to Y-values (upper and lower part of the scratched sheet)
%         ub=2*median(Y(sep==1))+64;lb=2*median(Y(sep==2))-max(Y)-64; % 
        nbins = 100;
        [n_cells_bin , center_bins] = hist(Y , nbins); % Number of cells in each bin
        cutoff_n_cells = 5; % Cutoff for the number of cells per bin
        lo_bin = find(n_cells_bin < cutoff_n_cells , 1 , 'first'); % First bin that has fewer than X cells
        hi_bin = find(n_cells_bin(lo_bin : end) > cutoff_n_cells , 1 , 'first') + lo_bin - 1; % First bin on the other side of scratch that has more than X cells
        lo_bound = mean(center_bins(lo_bin : lo_bin + 1)) + 20;
        hi_bound = mean(center_bins(hi_bin - 1 : hi_bin)) - 20;
        %{
        figure; hold on;
        hist(Y , nbins);
        plot([lo_bound , lo_bound] , [0 max(n_cells_bin)] , 'r--');
        plot([hi_bound , hi_bound] , [0 max(n_cells_bin)] , 'g--');
        %}        
        ub = lo_bound;
        lb = hi_bound;
        lo_up=lo(Y<ub & Y>(ub-strip)); % upper half of cells
        lo_dn=lo(Y>lb & Y<(lb+strip)); % lower half of cells
        lo_all=[lo_up; lo_dn]; % combines cells in the two strips above and beneath the scratch
        
        %lo_all=lo; %DG set to lo to ignore sorting of cells above and below scratch since not needed at 20x***CHANGE BACK
        i=0; % index for averaging between frames     
        numcell=length(lo_all);
        table(n).numbercells=numcell;
        
        % Saves the upper and lower bound of cell sheet
        table_bounds(n).filename = tempfile;
        lo_bound_matrix = NaN(numel(st:tstep:ed) , 1);
        hi_bound_matrix = lo_bound_matrix;
        for frame=st:tstep:ed
            i=i+1;
            tm1=[tracedata(lo_all,frame,1) tracedata(lo_all,frame,2)]; % x-y coordinates at time point t
            tm4=[tracedata(lo_all,frame+tstep,1) tracedata(lo_all,frame+tstep,2)]; % x-y coordinates at time point t+tstep
            d=tm4-tm1; % Vectors connecting two successive points. 
            [atot vtot]=cart2pol(d(:,1),d(:,2)); % polar coordinate transformation
            vtot=vtot(~isnan(vtot)); %remove NaNs
            vtot=vtot*pix*6/tstep; % scaling for IXµ with 6 frames/h (every 10min)
            v_ave(i)=mean(vtot);
            
            % Calculate bound of cells for each frame
            Y=tracedata(lo,frame,2);% Y-values of cells at ed.  %X=tracedata(lo,st,1);
            sep=kmeans(Y,2,'start',[1;max(Y)]); % seperates all cells in two groups according to Y-values (upper and lower part of the scratched sheet)
            %         ub=2*median(Y(sep==1))+64;lb=2*median(Y(sep==2))-max(Y)-64; %
            nbins = 100;
            [n_cells_bin , center_bins] = hist(Y , nbins); % Number of cells in each bin
            cutoff_n_cells = 5; % Cutoff for the number of cells per bin
            lo_bin = find(n_cells_bin < cutoff_n_cells , 1 , 'first'); % First bin that has fewer than X cells
            hi_bin = find(n_cells_bin(lo_bin : end) > cutoff_n_cells , 1 , 'first') + lo_bin - 1; % First bin on the other side of scratch that has more than X cells
            lo_bound = mean(center_bins(lo_bin : lo_bin + 1));
            hi_bound = mean(center_bins(hi_bin - 1 : hi_bin));
            lo_bound_matrix(i) = lo_bound;
            hi_bound_matrix(i) = hi_bound;
                 
        end %frame
        avgWellVel=mean(v_ave);
        stdWellDev=std2(v_ave);
        table(n).averageVelocity=avgWellVel;
        table(n).stdev=stdWellDev;
        
        table_bounds(n).lo = lo_bound_matrix;
        table_bounds(n).hi = hi_bound_matrix;
        %{
        figure;
        hold on;
        plot(lo_bound_matrix , 'r-');
        plot(hi_bound_matrix , 'b-');
        %}
        
        
        disp(tempfile); disp('is done onto next file!');
end
xtable = struct2table(table);
filewrite=[datadir, '\', date, 'CellSpeeds.xls'];
writetable(xtable,filewrite);
toc
disp('ALL DONE!');

save('table_bounds.mat' , 'table_bounds');
%% Histograms of first timepoints with edges saved as tifs
n=1;
  for n=1:length(filenames)
        tic
        tempfile=filenames{n,1}; %load filename
        table(n).filename=tempfile;
        load([rawdir tempfile]); %loads data from file
        lo=find(~isnan(tracedata(:,st,1))&~isnan(tracedata(:,ed,1))); % finds tracks going from st to ed
        
        Y=tracedata(lo,st,2);% Y-values of cells at ed.  %X=tracedata(lo,st,1);
        sep=kmeans(Y,2,'start',[1;max(Y)]); % seperates all cells in two groups according to Y-values (upper and lower part of the scratched sheet)
%         ub=2*median(Y(sep==1))+64;lb=2*median(Y(sep==2))-max(Y)-64; % 
        nbins = 100;
        [n_cells_bin , center_bins] = hist(Y , nbins); % Number of cells in each bin
        cutoff_n_cells = 5; % Cutoff for the number of cells per bin
        lo_bin = find(n_cells_bin < cutoff_n_cells , 1 , 'first'); % First bin that has fewer than X cells
        hi_bin = find(n_cells_bin(lo_bin : end) > cutoff_n_cells , 1 , 'first') + lo_bin - 1; % First bin on the other side of scratch that has more than X cells
        lo_bound = mean(center_bins(lo_bin : lo_bin + 1)) + 20;
        hi_bound = mean(center_bins(hi_bin - 1 : hi_bin)) - 20;
            temphist=figure; hold on;
            hist(Y , nbins); %make histogram
            plot([lo_bound , lo_bound] , [0 max(n_cells_bin)] , 'r--');
            plot([hi_bound , hi_bound] , [0 max(n_cells_bin)] , 'g--');
            saveas(temphist,[datadir,filesep,tempfile,'-histogram.tif']);
            hold off; close all;
  end
%% Plot
ID_plot=[2 7 6 5 4 3]; % [0 0.12 0.36 1.1 3.3]µM Y27632
scatter(1:6,v_mean_plate(2,ID_plot),'MarkerEdgeColor','b','MarkerFaceColor','b');hold on
scatter(1:6,v_mean_plate(3,ID_plot),'MarkerEdgeColor','r','MarkerFaceColor','r');
plot(1:6,v_mean_plate(2,ID_plot),'b');
plot(1:6,v_mean_plate(3,ID_plot),'r');hold off
axis([0 7 0 20]);

% %% Plot all healing traces
% for row=2:7
%     for col=1:12
%         y(:,1)=healing(row,col,:);
%         plot(1:72,y);hold on;
%     end
% end