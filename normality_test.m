% normality test
% intput: err = vecAll - vecAll_fit
H_list = [];
pValue_list = [];
mean_list = [];
ratio_list = [];
for d=1:2
    for x=1:5:500
        for y=1:5:500
            err1 = squeeze(err(d,x,y,1,:));
            err1 = interp1q((1:length(err1))', err1, (1:.1:length(err1))');
            if 1
                [H, pValue] = swtest(err1);
            else
                [H, pValue] = kstest(err1);
            end
            
            H_list(end+1) = H;
            pValue_list(end+1) = pValue;
            mean_list(end+1) = mean(err1);
            ratio_list(end+1) = mean(err1)/std(err1);
            
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
            
            