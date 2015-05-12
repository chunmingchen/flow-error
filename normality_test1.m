function [fit_ratio_ks, fit_ratio_chi, conf_95_ratio, pfield] = normality_test1(err, dof)

if nargin==1
    dof=1;
end

% normality test
% intput: err = vecAll - vecAll_fit
H_list1 = [];
pValue_list1 = [];
H_list2 = [];
pValue_list2 = [];

div = size(err,5)-dof;
conf_95_list = [];

pfield = 2*ones(size(err,1), size(err,2), size(err,3), size(err,4));

for z=1:size(err,4)
    for y=1:size(err,3)
        for d=1:size(err,1)
            for x=1:size(err,2)
                err1 = squeeze(err(d,x,y,z,:));
                s = sqrt(sum(err1.^2)/div); % standard error
                
                conf_95_list(end+1) = mean(abs(err1) < s*1.96);
               
                if s<1e-3
%                     H_list1(end+1) = 0;                    
                    continue
                end
                %pd = makedist('Normal', 'sigma', sqrt(mean(err1.^2)));
                pd = makedist('Normal', 'mu', 0, 'sigma', s);
                [H1, pValue1] = kstest(err1, 'CDF', pd);    
                [H2, pValue2] = chi2gof(err1, 'CDF', pd);
%                 [H2, pValue2] = chi2gof(err1);

hist(err1, 20)
                if ~isnan(pValue1)

                    H_list1(end+1) = H1;
                    pValue_list1(end+1) = pValue1;
                end
                if ~isnan(pValue2)

                    H_list2(end+1) = H2;
                    pValue_list2(end+1) = pValue2;
                end
                pfield(d,x,y,z) = pValue1;
                
%                 hist(err1);
%                 title(sprintf('chitest pvalue=%f, h=%f.  kstest pvalue=%f, h=%f', pValue2, H2, pValue1, H1))
%                 pause
                
            end
        end
    end
end
% hist(pValue_list,20);
disp(sprintf('length: ks=%d chi=%d', length(H_list1), length(H_list2)))
fit_ratio_ks=1-sum(H_list1)/length(H_list1)     ;
fit_ratio_chi=1-sum(H_list2)/length(H_list2)     ;
conf_95_ratio = sum(conf_95_list)/length(conf_95_list)
hist(conf_95_list,20)
end