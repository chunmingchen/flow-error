function [fit_ratio_ks, fit_ratio_chi] = normality_test0(err)

% normality test
% intput: err = vecAll - vecAll_fit
H_list1 = [];
pValue_list1 = [];
H_list2 = [];
pValue_list2 = [];
for z=1:size(err,4)
    for y=1:size(err,3)
        for d=1:size(err,1)
            for x=1:size(err,2)
                err1 = squeeze(err(d,x,y,1,:));
                %pd = makedist('Normal', 'mu', mean(err1), 'sigma', std(err1));
                %[H1, pValue1] = kstest(err1, 'CDF', pd);    
                %[H2, pValue2] = chi2gof(err1, 'CDF', pd);
                [H1, pValue1] = kstest(err1);    
                [H2, pValue2] = chi2gof(err1);

                if ~isnan(pValue1)

                    H_list1(end+1) = H1;
                    pValue_list1(end+1) = pValue1;
                end
                if ~isnan(pValue2)

                    H_list2(end+1) = H2;
                    pValue_list2(end+1) = pValue2;
                end
            end
        end
    end
end
% hist(pValue_list,20);
disp(sprintf('length: ks=%d chi=%d', length(H_list1), length(H_list2)))
fit_ratio_ks=1-sum(H_list1)/length(H_list1)     ;
fit_ratio_chi=1-sum(H_list2)/length(H_list2)     ;

end