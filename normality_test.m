% normality test
% intput: err = vecAll - vecAll_fit
H_list = [];
pValue_list = [];
mean_list = [];
ratio_list = [];
for d=1:2
    for x=1:10:data_w
        for y=1:10:data_h
            err1 = squeeze(err(d,x,y,1,:));
%             err1 = interp1q((1:length(err1))', err1, (1:.1:length(err1))');
            if 0
                [H, pValue] = swtest(err1);
            elseif 0
                [H, pValue] = kstest(err1);
            elseif 0
%                 criticalvalue = sum(err
            else
                pd = makedist('Normal', 'sigma', sqrt(mean(err1.^2)));
%                 pd = fitdist(err1, 'Normal');
                if 1
                    [H, pValue] = kstest(err1, 'CDF', pd);    
                elseif 0
                    [H, pValue] = adtest(err1, 'CDF', pd);
                else
                    [H, pValue] = chi2gof(err1, 'CDF', pd);
                end
            end
            
            if ~isnan(pValue)
            
                H_list(end+1) = H;
                pValue_list(end+1) = pValue;
                mean_list(end+1) = mean(err1);
                ratio_list(end+1) = mean(err1)/std(err1);
            end
%             if (H==1)
%                 hist(err1, 20)
%                 pValue
%                 x
%                 y
%                 d
%                 
%                 pause
%             end
%                 normplot(err1)
%                 pause
        end
    end
end
hist(pValue_list,20);
fit_ratio=1-sum(H_list)/length(H_list)     
            