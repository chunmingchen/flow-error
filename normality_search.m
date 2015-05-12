function [ks_bezier_search, chi_bezier_search, conf_ratio_bezier_search, pfield_bezier_search, err_rms_bezier_search] = normality_search(vec)
PLOT=1
% normality test
% intput: err = vecAll - vecAll_fit
H_list1 = [];
pValue_list1 = [];
H_list2 = [];
pValue_list2 = [];

conf_95_list = [];

pfield_bezier_search = 2*ones(size(vec,1), size(vec,2), size(vec,3), size(vec,4));
err_rms_bezier_search= [];
cases=0;
for z=1:size(vec,4)
    z
    for y=1:size(vec,3)
        for d=1:size(vec,1)
            for x=1:size(vec,2)
                cases = cases+1;
                vec1 = squeeze(vec(d,x,y,z,:))';
                N=length(vec1);
                div = N-1;
                tt = 0:(N-1);
                [yy0, ctrl0] = bezierfit1(tt, vec1, tt);
                TESTS=31;
                plist1=zeros(TESTS,1);
                plist2=zeros(TESTS,1);
                count=1;
                ctrls = linspace(ctrl0-5, ctrl0+5, TESTS); % tests
                if PLOT
                    subplot(2,1,1)
                    plot(tt, vec1, 'xg')
                    hold on
                end
                for ctrl = ctrls
                    yy =  ((N-1-tt)/(N-1)).^2 * vec1(1) + 2 * (N-1-tt).*tt/(N-1)^2 .* ctrl + tt.^2/(N-1)^2*vec1(N);

                    err1 = yy-vec1;
                    if PLOT
                        plot(tt, yy)
                    end
                    s = sqrt(sum(err1.^2)/div); % standard error
                    pd = makedist('Normal', 'mu', 0, 'sigma', s);
                    [H1, pValue1] = kstest(err1, 'CDF', pd);    
                    [H2, pValue2] = chi2gof(err1, 'CDF', pd);
    %                 [H2, pValue2] = chi2gof(err1);
    
                    plist1(count)=pValue1;
                    plist2(count)=pValue2;
                    count = count+1;
                end
                if PLOT
                    plot(tt, vec1, 'xg-')                
                    hold off
                    subplot(2,1,2)
                    plot(ctrls-ctrl0, plist1, ctrls-ctrl0, plist2)
                    legend('KS', 'CHI2')
                    saveas(gca, sprintf('fig/bezierfit_case%d.png', cases))
%                     pa[use
                end
                
                % get max
                [pmax1, I] = max(plist1);
                ctrlmax = ctrls(I);
                pmax2 = plist2(I);
                yymax =  ((N-1-tt)/(N-1)).^2 * vec1(1) + 2 * (N-1-tt).*tt/(N-1)^2 .* ctrlmax + tt.^2/(N-1)^2*vec1(N);
                err1 = yymax-vec1;
                
                % statistics
                conf_95_list(end+1) = mean(abs(err1) < s*1.96);
               
                
                %hist(err1, 20)
                pValue_list1(end+1) = pmax1;
                H_list1(end+1) = pmax1 > 0.05;
                pValue_list2(end+1) = pmax2;
                H_list2(end+1) = pmax2 > 0.05;
                pfield_bezier_search(d,x,y,z) = pmax1;
    
                err_rms_bezier_search(end+1) = sqrt(mean(err1.^2));
    %                 hist(err1);
    %                 title(sprintf('chitest pvalue=%f, h=%f.  kstest pvalue=%f, h=%f', pValue2, H2, pValue1, H1))
    %                 pause
            end
        end
    end
end
% hist(pValue_list,20);
disp(sprintf('length: ks=%d chi=%d', length(H_list1), length(H_list2)))
ks_bezier_search=1-sum(H_list1)/length(H_list1)     ;
chi_bezier_search=1-sum(H_list2)/length(H_list2)     ;
conf_ratio_bezier_search = sum(conf_95_list)/length(conf_95_list)
%hist(conf_95_list,20)
end