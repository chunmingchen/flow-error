n=1024

x=0:20/(n-1):20;

% y=sin(x);
y=sin(x).*cos(x*2);
% y=sin(x).*cos(x*2).*x + randn(size(x))*.01;
% y=randn(size(x));
% y=rand(size(x));

t_ary = (1:n);

if 0 % bezier
%     v = bezierfit1((t_ary-1), y, (t_ary-1))   ;
    v = bezierfitquartic((t_ary-1), y, (t_ary-1))' ;
%     v = bezierfitQuintic((t_ary-1), y, (t_ary-1))' ;
    subplot(2,2,1);
    plot(t_ary, y, t_ary, v);
elseif 0
    p = polyfit(t_ary, y, 24); 
    v = polyval(p, t_ary);

elseif 1
    v = myhaar(y', 9)';
    
else
    [af, sf] = farras; 
    w = dwt(y,1,af); 
    v = idwt(w,1,sf); 
    
end

subplot(2,2,1);
plot(t_ary, y, t_ary, v);
title('Signal')

err = y-v;
% err = y;
subplot(2,2,2)
plot(err, '+')
title ('Error')

s = sqrt(sum(err.^2)/(n-1)) % standard error

subplot(2,2,3)
hist(err,25)
hold on
% plot(-1:.1:1, gauss_distribution(-1:.1:1, 0, s) * n/10, 'r')
hold off


conf_95 = mean(abs(err) < s*1.96)
pd = makedist('Normal', 'mu', 0, 'sigma', s);
[H1, pValue1] = kstest(err, 'CDF', pd);    
disp('kstest')
[H1, pValue1]
[H1, pValue1] = kstest(err);    
[H1, pValue1]
[H2, pValue2] = chi2gof(err, 'CDF', pd);
disp('chi2')
[H2, pValue2]
[H3, pValue3] = chi2gof(err);
disp('chi2-natural')
[H3, pValue3]
title (sprintf('Normality: KS: %f, Chi2: %f', pValue1, pValue2))

subplot(2,2,4)
qqplot(err)


% exponential
pd = fitdist(abs(err)','gamma')
[H1, pValue1] = kstest(abs(err), 'CDF', pd);    
disp('exponential kstest')
[H1, pValue1]

[count, centers] = hist(abs(err), 25);
subplot(2,2,3)
hist(abs(err), 50)
hold on 
y = pdf(pd,centers);
%plot(centers,y*n/10,'r')
hold off
subplot(2,2,4)
qqplot(abs(err), pd)


