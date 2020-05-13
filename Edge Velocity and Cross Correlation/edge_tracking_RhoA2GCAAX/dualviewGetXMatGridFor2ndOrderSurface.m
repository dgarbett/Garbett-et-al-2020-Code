function xMat=dualviewGetXMatGridFor2ndOrderSurface(xvals,yvals)

[x,y]=meshgrid(xvals,yvals);
x=x(:); y=y(:);
xMat=zeros(length(x(:)),5);
xMat(:,1)=x.*x;
xMat(:,2)=y.*y;
xMat(:,3)=x.*y;
xMat(:,4)=x;
xMat(:,5)=y;
