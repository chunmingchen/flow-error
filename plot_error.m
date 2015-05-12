xvalues = { [12 24 47], [12, 24], [25 50 100] };
yvalues = {    [3.1381	4.2965	6.9277
2.2047	2.9946	4.8209
1.9845	2.6719	4.3978], 
[0.3859	0.44689
0.329	0.37171
0.3055	0.33798],
[3.931	4.4094	4.921197507
2.9699	3.6735	4.468828321
2.5788	3.2462	3.952123677],
[-0.0147	-0.3792	-0.2658
0.0066	-0.0974	-0.2759
0.0034	0.0134	-0.0239],
[6.05E-03	-0.00032762
7.22E-03	0.0022336
-9.68E-04	-0.00063945],
[-2.92E-01	-0.052825	-0.15951
-1.24E-01	-0.21459	-0.3178
-1.74E-02	-0.037677	-0.054516]
};
%legends = {{'Linear Interpolation', 'Linear Interp. w/ Dbl. Sampling', 'Quadratic Bezier Fit'}};
legends = {{'LI', 'LI-DS', 'QB'}};
xlabels = {'Sampling Interval (Time Steps)'};
ylabels = {'RMSE', 'RMSE', 'RMSE', 'Mean Error', 'Mean Error', 'Mean Error'};
titles = {'Isabel', 'Plume', 'Climate', 'Isabel', 'Plume', 'Climate'};
base_path = '~/Dropbox/0osuflow/Paper/LDAV14/fig/matlab';


close all
jet_inv = jet;
jet_inv = jet_inv(length(jet_inv):-1:1,:)
% RMSE
for i=1:length(xvalues)
    figure;
    colormap(jet_inv)
    
    yvalues{i} = abs(yvalues{i});
    bar(yvalues{i}');
    if i==1
        h = legend(legends{1}, 'Location', 'NorthWest');
        set(h, 'FontSize', 28.0);
        set(h, 'box', 'off');
    end
    set(gca, 'FontSize', 36.0);
    xlabel(xlabels{1})
    ylabel(ylabels{i})
    title(titles{i})
    xlim([0.5,size(yvalues{i}, 2)+.5])
    ylim([0, max(max(yvalues{i}))*1.1])
    set(gca, 'XTickMode', 'manual')
    set(gca,'XTick', 1:(xvalues{i}))
    set(gca,'XTickLabel',xvalues{i})
    saveas(gca, sprintf('%s/%s_%s.eps', base_path, ylabels{i}, titles{i}), 'psc2');
end


% Error
for i=1:length(xvalues)
    figure;
%    yvalues{i} = abs(yvalues{i});
    hold on
    bar(yvalues{i}');
    bar((yvalues{i+3}'), 'w');
    hold off
    %h = legend(legends{1});
    %set(h, 'FontSize', 12.0);
    set(gca, 'FontSize', 36.0);
    xlabel(xlabels{1})
    ylabel('Error')
    title(titles{i})
    xlim([0.5,size(yvalues{i}, 2)+.5])
    ylim([min(0, min(min(yvalues{i+3}))*1.1), max(max(yvalues{i}))*1.1])
    set(gca, 'XTickMode', 'manual')
    set(gca,'XTick', 1:(xvalues{i}))
    set(gca,'XTickLabel',xvalues{i})
    saveas(gca, sprintf('%s/Error_%s.eps', base_path,  titles{i}), 'psc2');
end

% Mean Error  percentage
for i=1:length(xvalues)
    figure;
%    yvalues{i} = abs(yvalues{i});
    bar(100*abs(yvalues{i+3}(3,:) ./ yvalues{i}(3,:)));
    if i==1
        h = legend(legends{1}{3}, 'Location', 'NorthWest');
        set(h, 'FontSize', 28.0);
        set(h, 'box', 'off')
    end
    set(gca, 'FontSize', 36.0);
    xlabel(xlabels{1})
    h = ylabel('Mean Error / Standard Error (%)')
    set(h, 'FontSize', 27.0);
    title(titles{i})
    xlim([0.1,size(yvalues{i}, 2)+.9])
    ylim([0 1.5])
    set(gca, 'XTickMode', 'manual')
    set(gca,'XTick', 1:(xvalues{i}))
    set(gca,'XTickLabel',xvalues{i})
    
    saveas(gca, sprintf('%s/MeanErrorP_%s.eps', base_path,  titles{i}), 'psc2');
end
