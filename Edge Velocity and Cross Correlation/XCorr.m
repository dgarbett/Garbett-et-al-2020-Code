%%
files = dir('*4color*.mat'); %set to output from edge velocity analysis
frontPercent = 5;
%this is correct code
window = [42:57];
%signal1 = squeeze(protvalsWindowF(window,:));
%signal2 = squeeze(edgeSensorVals(:,1:end-1,));
    
corrData = [];

shifts = -30:30;
indexes = [];
%for window = 1:80
  for i=1:length(shifts)
      numFrames = 298;
      indF1=max(1,1+ceil(shifts(i))); indF2=min(numFrames,numFrames+ceil(shifts(i)));
      indP1=max(1,ceil(-1*shifts(i))); indP2=min(numFrames-1,numFrames+floor(-1*shifts(i)));
      
      %signal1 = squeeze(protvalsWindowF(window,:));
      signal1 = squeeze(edgeSensorVals(window,1:end-1,1));
      for s = 1:size(edgeSensorVals,3)
          signal2 = squeeze(edgeSensorVals(window,1:end-1,s));
          if shifts(i) == 0
              [temp, ptemp] = myNanCorrcoef(vect(signal1),vect(signal2));
              %temp = corrc(1,2);
              indexesF1(i,1:2) = [1, size(signal1,2)];
              indexesP1(i,1:2) = [1, size(signal2,2)];
          else
              
              [temp, ptemp]=myNanCorrcoef(vect(signal1(:,indF1:indF2)),vect(signal2(:,indP1:indP2)));
              indexesF1(i,1:2) = [indF1, indF2];
              indexesP1(i,1:2) = [indP1, indP2];
          end
         % corrData(window,i,s) = temp;
         corrData(i,s) = temp;
      end
  end
%end;

%%
figure;
yyaxis left
plot(shifts*1/12,corrData(:,3), 'k-o','DisplayName','TP/FT');
ylabel('Cross Correlation MPAct')
hold on;
yyaxis right
plot(shifts*1/12,corrData(:,1), 'r-o','DisplayName','FT');
legend;
xlabel('Lag relative to protrusion (min)');
ylabel('Cross Correlation prot')
%plot(shifts*1/2,corrData(:,4), 'b-o','DisplayName','MPAct(d1-6)');

%%
 %show the shifts graphically
  
yF1 = -30:30;
yP1 = yF1+.25;
figure;
hold on;
for i = 1:length(shifts)
    plot(indexesF1(i,:),[yF1(i) yF1(i)], 'b-','LineWidth',4);
    plot(indexesP1(i,:),[yP1(i) yP1(i)], 'r-','LineWidth',4);
end

