% compute e = X - f(X)
% X is distributed as (xi, bincount)
function val = sample_histogram(samples, edge, bincount)

assert(length(edge)==length(bincount));
% pdf = bincount/sum(bincount);
% cdf = cumsum(pdf);

r=length(bincount):-1:1;
revbincount = bincount(r);
revpdf = revbincount/sum(revbincount);
revcdf = cumsum(revpdf);
cdf = 1-revcdf(r);

uniform_x = rand(samples,1);
val = zeros(samples,1);
for i=1:samples
    bin = max(find(uniform_x(i) > cdf));
    if isempty(bin)
        bin=1;
    end
    % bin: edge(bin) ~ edge(bin+1)
    if (bin<length(edge))
        val(i) = edge(bin) + (edge(bin+1)-edge(bin))*rand(1);
    else
        val(i)=edge(end);
    end
    
end

% assert(size( edge,1)==size(bincount,1) && size( edge,2)==1 )
% 
% r=length(bincount):-1:1;
% 
% revbincount = bincount(r);
% revpdf = revbincount/sum(revbincount);
% revcdf = cumsum(revpdf);
% cdf = 1-revcdf(r);
% uniform_x = rand(samples,1);
% val = interp1q([0;cdf], [ edge(1); edge], uniform_x);


end