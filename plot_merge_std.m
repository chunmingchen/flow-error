xvalues = { [12 24 47], [12, 24], [25 50 100] };
yvalues = {   [0.8393	0.6397	0.8397
0.9911	0.9028	0.8657
0.9992	0.9669	0.8878]
[0.6518	0.4657
0.9112	0.7587
0.984	0.9669]
[0.8092	0.4449	0.4388
0.9547	0.6172	0.7907
0.9944	0.8397	0.9302]
};
legends = {{'1x STD', '2x STD', '3x STD'}};
xlabels = {'Sampling Interval (Time Steps)'};
ylabels = {'Coverage (%)'};
titles = {'Isabel', 'Plume', 'Climate'};
base_path = '~/Dropbox/0osuflow/Paper/LDAV14/fig/matlab';


close all
jet_inv = jet;
jet_inv = jet_inv(length(jet_inv):-1:1,:)
% RMSE
for i=1:length(xvalues)
    figure;
    colormap (winter)
    
    bar(100*yvalues{i}');
    if i==1
         h = legend(legends{1}, 'Location', 'SouthEast');
         set(h, 'FontSize', 24.0);
%           set(h,'box','off') 
    end
    set(gca, 'FontSize', 36.0);
    xlabel(xlabels{1})
    ylabel(ylabels{1})
    title(titles{i})
    xlim([0.5,size(yvalues{i}, 2)+.5])
%     ylim([0, 100])
%     set(gca,'YTick', 0:10:100);
%     set(gca,'YTickLabel', {0, '', 20, '', 40, '', 60, '', 80, '', 100});

%     if i==3
        set(gca, 'YTick', 0:20:100);
%     end
    
    ylim([0, 101])
    set(gca, 'XTickMode', 'manual')
    set(gca,'XTick', 1:length(xvalues{i}))
    set(gca,'XTickLabel',xvalues{i})
%     my_xticklabels(gca, 1:length(xvalues{i}), xvalues{i})
    saveas(gca, sprintf('%s/%s_%s.eps', base_path, 'std', titles{i}), 'psc2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

yvalues = {   [0.99	0.8947	0.8637]
    [0.9038	0.749]
       [0.9505	0.6086	0.7776]
};
legends = {{'1.96x STD'}};



% RMSE
for i=1:length(xvalues)
    figure;
    colormap (jet)
    
    bar(100*yvalues{i}');
    if i==1
         h = legend(legends{1}, 'Location', 'SouthEast');
         set(h, 'FontSize', 24.0);
%           set(h,'box','off') 
    end
    set(gca, 'FontSize', 36.0);
    xlabel(xlabels{1})
    ylabel(ylabels{1})
    title(titles{i})
    xlim([0.5,size(yvalues{i}, 2)+.5])
%     ylim([0, 100])
%     set(gca,'YTick', 0:10:100);
%     set(gca,'YTickLabel', {0, '', 20, '', 40, '', 60, '', 80, '', 100});

%     if i==3
        set(gca, 'YTick', 0:20:100);
%     end
    
    ylim([0, 101])
    set(gca, 'XTickMode', 'manual')
    set(gca,'XTick', 1:length(xvalues{i}))
    set(gca,'XTickLabel',xvalues{i})
%     my_xticklabels(gca, 1:length(xvalues{i}), xvalues{i})
    saveas(gca, sprintf('%s/%s_%s.eps', base_path, 'std_1.96', titles{i}), 'psc2');
end
