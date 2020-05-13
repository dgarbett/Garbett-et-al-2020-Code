function imnew=CropJitter(imold,cx1,cx2,cy1,cy2,dx,dy)

% cx1,cx2,cy1,cy2: margin
% dx,dy: jitter

%% get floor and remain
dxf=floor(dx);dxr=dx-dxf;
dyf=floor(dy);dyr=dy-dyf;

%% get crude image
new11=imold(1+cy1-dyf:end-cy2-dyf,1+cx1-dxf:end-cx2-dxf);
new12=imold(cy1-dyf:end-cy2-dyf-1,1+cx1-dxf:end-cx2-dxf);
new21=imold(1+cy1-dyf:end-cy2-dyf,cx1-dxf:end-cx2-dxf-1);
new22=imold(cy1-dyf:end-cy2-dyf-1,cx1-dxf:end-cx2-dxf-1);

%% intrapolation
new1y=(1-dyr)*new11+dyr*new12;
new2y=(1-dyr)*new21+dyr*new22;
imnew=(1-dxr)*new1y+dxr*new2y;