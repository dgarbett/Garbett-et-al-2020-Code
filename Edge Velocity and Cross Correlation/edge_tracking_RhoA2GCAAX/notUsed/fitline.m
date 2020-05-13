function res = fitline(t,y)
%FITLINE returns the best fit line to the data, as well as an r^2 value

t=t(:);
y=y(:);
t3=t(~isnan(t));
y3=y(~isnan(t));
t2=t3(~isnan(y3));
y2=y3(~isnan(y3));
x = [t2 ones(size(t2))];
a = x\y2;
res.m = a(1);
res.b = a(2);
if length(y)>2
    temp = corrcoef(t2,y2);
    res.r = temp(1,2);
    res.r2 = res.r*res.r;
else
    res.r=NaN;
    res.r2=NaN;
end
