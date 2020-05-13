function protrusionVals=computeProtrusionValues(edgeCoors,edgeCoorsSmoothed,edgeCoors2)
%This function measures the local protrusion or retraction of the cell edge
%adjacent to each of the points specified in edgeCoors

% map edgeCoors to mask2 and get normal vectors
normV=computeNormalVectorsFromParametrizedCellEdge(edgeCoorsSmoothed);
% figure; hold on; plot(edgeCoors(:,1),edgeCoors(:,2)); quiver(edgeCoors(:,1),edgeCoors(:,2),0.2*normV(:,1),0.2*normV(:,2));

% measure protrusion for each point in edgeCoors (orthogonal to cell edge)
protrusionVectors=edgeCoors2-edgeCoors;
protrusionVals=sum(protrusionVectors.*normV,2); % dot product between protrusionVectors and normV, "projects" protrusionVectors onto normV.

% figure; hold on;
% plot(edgeCoors(:,1),edgeCoors(:,2),'.-','Color','b');
% plot(edgeCoors2(:,1),edgeCoors2(:,2),'.-','Color','k');
% % scatter(edgeCoors(protrusionVals>5,1),edgeCoors(protrusionVals>5,2),'r');
% % scatter(edgeCoors(protrusionVals<-5,1),edgeCoors(protrusionVals<-5,2),'b');
% for i=1:size(edgeCoors,1)
%     plot([edgeCoors(i,1) edgeCoors2(i,1)],[edgeCoors(i,2) edgeCoors2(i,2)],'.-');
% end