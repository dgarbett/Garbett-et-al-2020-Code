function d=distanceFromPointsToLine(points, line)
%returns the distance between each point and a line using the formula obtained
%from http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line (in
%two dimensions). Points should be a matrix where each row is one point and line should be
%a with the first two numbers indicating a point on the line, and the last
%two indicating a vector for the direction of the line.

a=line(1:2); a=a(:)';
n=line(3:4); n=n/norm(n); n=n(:)';

a=repmat(a,[size(points,1) 1]);
n=repmat(n,[size(points,1) 1]);
vectors= (a - points) - n.*repmat(sum((a-points).*n,2),[1 2]);

d=sqrt(vectors(:,1).*vectors(:,1) + vectors(:,2).*vectors(:,2));
