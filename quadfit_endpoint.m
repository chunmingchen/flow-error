% fit column vectors
function [errstd errsum yy a b c rmserr, err] = quadfit_endpoint(x,y, xx)

sum_y  = sum(y);
sum_y2 = sum(y.*y);
sum_xy = sum(x.*y);
sum_x2y= sum(x.*x.*y);
sum_x  = sum(x);
sum_x2 = sum(x.*x);
sum_x3 = sum(x.^3);
sum_x4 = sum(x.^4);

if length(x)~=length(y)
    error('length not equal');
end
n = length(x);
yn = y(end);
y1 = y(1);
xn = x(end);
x1 = x(1);
if (x1 ~= 0)
    error ('x1 must be zero')
end

a=(sum_x3*(yn-y1)/xn - sum_x2y - yn*sum_x2 + 2*y1*sum_x2 + xn*sum_xy - xn*y1*sum_x) ...
							/ (2*sum_x3*xn-sum_x4-sum_x2*xn*xn);
b= ((yn-y1)-a*(xn*xn))/xn;
c= y(1);

yy = a*(x.*x)+b*x+c;
err = yy-y;
errsum = abs(sum(err))
errstd = std(err)
rmserr = sqrt(mean(err.*err));

if nargin>=3 
    yy =  a*(xx.*xx)+b*xx+c;
end

end