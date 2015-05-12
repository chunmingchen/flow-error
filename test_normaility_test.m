pvalue1_list=[];
pvalue2_list=[];
err1=randn(1,1000);
for i=1:1
    mean(err1)
    std(err1)
    pd = makedist('Normal', 'mu', 0, 'sigma', 1);
    [H1, pValue1] = kstest(err1, 'CDF', pd)
    [H2, pValue2] = chi2gof(err1, 'CDF', pd)
    % [H2, pValue2] = chi2gof(err1)
    pvalue1_list(end+1)=pValue1;
    pvalue2_list(end+1)=pValue2;
    err1 = [err1 err1];
    hist(err1,30)
    set(gca, 'YTick', []);
    set(gca, 'FontSize', 30.0);
    hold off
    x=-4:.1:4;
    plot(x,pdf(pd,x)*1000,'r','LineWidth',2)
    hold off
%     pause
    saveas(gca, 'dist2.png')
end
% hist(pvalue1_list,100)