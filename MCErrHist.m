% compute e = X - f(X)
% X is distributed as (xi, bincount)
function estErrDist = MCErrHist(xi, bincount, quadfun, err_xi)
    
assert(length(xi)==length(bincount));
pdf = bincount/sum(bincount);
cdf = cumsum(pdf);

SAMPLES = 10000;
uniform_x = rand(SAMPLES,1);
dist_x = zeros(SAMPLES,1);
for i=1:SAMPLES
    bin = max(find(uniform_x(i) > cdf));
    if isempty(bin)
        dist_x(i) = xi(1);
    else
        dist_x(i) = xi(bin);
    end
end
f_dist_x = quadfun(1).*dist_x.*dist_x + quadfun(2).*dist_x + quadfun(3);
estErrBincount = histc(dist_x - f_dist_x , err_xi);
estErrDist = estErrBincount / SAMPLES;


end