% fit column vectors
% [x1 x2 x3 x4 ...]
% [y1 y2 y3 y4 ...]
function [errstd errsum yy ctrl err] = bezierfit(t,y, tt)

if t(1)~=0
    error('t(1)~=0')
end

N = length(t);
if (size(y,2)~=N)
    error ('size(y,2)~=N')
end
ctrl=zeros(size(y,1),1);
for d=1:size(y,1)
    ctrl(d) = sum( (N-1-t).*t.*y(d,:)*(N-1)^2 ) - sum( (N-1-t).^3 .* t .* y(d,1)) - sum( (N-1-t).*t.^3.*y(d,N) ) ;
    ctrl(d) = ctrl(d) / 2 / sum( (N-t-1).^2 .* t.^2 );
end

errstd=0
errsum=0

yy=zeros(size(y,1),length(tt));
if nargin>=3 
    for d=1:size(y,1)
        yy(d,:) =  ((N-1-tt)/(N-1)).^2 * y(d,1) + 2 * (N-1-tt).*tt/(N-1)^2 .* ctrl(d) + tt.^2/(N-1)^2*y(d,N);
    end
end

end