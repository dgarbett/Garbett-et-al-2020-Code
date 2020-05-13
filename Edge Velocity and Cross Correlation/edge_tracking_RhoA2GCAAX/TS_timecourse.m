clear;clc;close all; 
root = '/Volumes/BisackUp/Meyer_Lab/CDH5-TS-Thrombin';
datadirTS=[root,filesep,'CDH5-TS-1min_1',filesep,'data'];
datadirTL=[root,filesep,'CDH5-TSTL-1min_1',filesep,'data'];

%%

timepoints = 25;
filenamesTS=getFilenames(datadirTS);
filenamesTL=getFilenames(datadirTL);

filenamesTS=filenamesTS(boolRegExp(filenamesTS,'RatioData.mat'));
filenamesTL=filenamesTL(boolRegExp(filenamesTL,'RatioData.mat'));
medIntensity = cell(timepoints,1);
    
for f = 1:length(filenamesTS)
    load([datadirTS,filesep,filenamesTS{1}]);
    medIntensityTemp = cell(length(maskFinal),1);
    for t = 1:length(maskFinal)
        CC = bwconncomp(maskFinal{t});
        stats = regionprops(CC,imRatio{t},'PixelValues');
        medIntensityTemp{t,1} = [medIntensity{t,1} cellfun(@(x) prctile(x,5),struct2cell(stats))];
    end
    
    %look at average value before thrombin;
    medIntensitySum = cellfun(@sum,medIntensityTemp);
    medIntensityNum = cellfun(@length,medIntensityTemp);
    averagePreThrom =  sum(medIntensitySum(1:10))/sum(medIntensityNum(1:10));
    
    
    for t=1:length(maskFinal)
         medIntensity{t,1} = [medIntensity{t,1} medIntensityTemp{t,1}/averagePreThrom];
    end;
end

medVec = cell2mat(medIntensity')';
groupVec =[];
for t = 1:length(medIntensity);
    groupVec = [groupVec ;repmat(t,length(medIntensity{t,1}),1)];
end

boxplot(medVec,groupVec)
%title('CDH5-TS 
%scatter(groupVec,medVec,12,'jitter','on','jitterAmount',0.5);
%medIntensityTL = cell(timepoints,1);
    
for f = 1:length(filenamesTL)
    load([datadirTL,filesep,filenamesTL{1}]);
     medIntensityTemp = cell(length(maskFinal),1);
    for t = 1:length(maskFinal)
        CC = bwconncomp(maskFinal{t});
        stats = regionprops(CC,imRatio{t},'PixelValues');
        medIntensityTemp{t,1} = [medIntensity{t,1} cellfun(@(x) prctile(x,5),struct2cell(stats))];      
    end
    medIntensitySum = cellfun(@sum,medIntensityTemp);
    medIntensityNum = cellfun(@length,medIntensityTemp);
    averagePreThrom =  sum(medIntensitySum(1:10))/sum(medIntensityNum(1:10));
    
    for t=1:length(maskFinal)
         medIntensityTL{t,1} = [medIntensityTL{t,1} medIntensityTemp{t,1}/averagePreThrom];
    end;
    
end

medVecTL = cell2mat(medIntensityTL')';
groupVecTL =[];
for t = 1:length(medIntensityTL);
    groupVecTL = [groupVecTL ;repmat(t,length(medIntensityTL{t,1}),1)];
end
figure(2)
boxplot(medVecTL,groupVecTL)
%title('CDH5-TS 
%scatter(groupVecTL,medVecTL,12,'filled','jitter','on','jitterAmount',0.5);



