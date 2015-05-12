xvalues = { [12 24 47], [12, 24], [25 50 100] };
yvalues = {   [16.1315	36.1786	58.4566
9.0292	48.1381	81.0569
4.5454	12.7662	31.9445]
[0.5276	1.3891
1.1474	1.1248
0.1249	0.379]
[19.1771	38.764	71.7922
24.8951	63.5604	121.9786
11.2699	28.8575	60.8385]
};
%legends = {{'Directly Connect', 'RK4 Tracing Double Sampling', 'Forward-backward Intersecting       '}};
legends = {'Directly Connect', 'Fw-DS', 'Fw-Bw'};
xlabels = {'Sampling Interval (Time Steps)'};
ylabels = {'Average Error (Voxels)'};
titles = {'Isabel', 'Plume', 'Climate'};
base_path = '~/Dropbox/0osuflow/Paper/PacificVis15/fig/matlab';

% swap
tmp = legends{1}; legends{1}=legends{2}; legends{2}=tmp;
for i=1:3
    tmp = yvalues{i}(1,:); yvalues{i}(1,:)=yvalues{i}(2,:); yvalues{i}(2,:)=tmp;
end

close all
jet_inv = jet;
jet_inv = jet_inv(length(jet_inv):-1:1,:)
% RMSE
for i=1:length(xvalues)
    figure;
    colormap(jet_inv)
    
    bar(yvalues{i}');
    if i==1
         h = legend(legends, 'Location', 'NorthWest');
         set(h, 'FontSize', 28.0);
          set(h,'box','off') 
    end
    set(gca, 'FontSize', 36.0);
    xlabel(xlabels{1})
    ylabel(ylabels{1})
    title(titles{i})
    xlim([0.5,size(yvalues{i}, 2)+.5])
%     ylim([0, 100])
%     set(gca,'YTick', 0:10:100);
%     set(gca,'YTickLabel', {0, '', 20, '', 40, '', 60, '', 80, '', 100});

    if i==3
        set(gca, 'YTick', 0:50:200);
    end
    
    ylim([0, max(max(yvalues{i}))*1.1])
    set(gca, 'XTickMode', 'manual')
    set(gca,'XTick', 1:length(xvalues{i}))
    set(gca,'XTickLabel',xvalues{i})
%     my_xticklabels(gca, 1:length(xvalues{i}), xvalues{i})
    saveas(gca, sprintf('%s/%s_%s.eps', base_path, 'pathline_err', titles{i}), 'psc2');
end


