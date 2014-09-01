% xvalues = { {{'Chi2';'12'}, {'KS';'12'}, {'Chi2';'24'} ,{'KS';'24'}, {'Chi2';'47'}, {'KS';'47'} },  ...
%     {{'Chi2';'12'} {'KS';'12'} {'Chi2';'24'} , {'KS';'24'}}, ...
%     { {'Chi2';'25'} {'KS';'25'} {'Chi2';'50'} {'KS';'50'} {'Chi2';'100'}  {'KS';'100'}} };
xvalues = { {{'KS';'12'}, {'Chi2';'12'}, {'KS';'24'} ,{'Chi2';'24'}, {'KS';'47'}, {'Chi2';'47'} },  ...
    {{'KS';'12'} {'Chi2';'12'} {'KS';'24'} , {'Chi2';'24'}}, ...
    { {'KS';'25'} {'Chi2';'25'} {'KS';'50'} {'Chi2';'50'} {'KS';'100'}  {'Chi2';'100'}} };
yvalues = {   [0.5256	0.6264	0.308	0.3736	0.1088	0.1648
0.752	0.8344	0.528	0.6144	0.2704	0.3568
1	0.9976	0.9896	0.9736	0.8848	0.84]
[0.23333	0.67333	0.30667	0.36333
0.92333	0.98	0.44	0.63
1	1	1	0.99]
[0.33816	0.37867	0.25965	0.31848	0.2157	0.2251
0.50559	0.57446	0.38368	0.43094	0.2508	0.2683
0.99306	0.9755	0.93268	0.92052	0.7942	0.7427]
};
legends = {{'Linear Interpolation', 'Linear Interp. w/ Double Sampling', 'Quadratic Bezier Spline'}};
xlabels = {'Test Method and Sampling Interval'};
ylabels = {'Rejection Rate (%)'};
titles = {'Isabel', 'Plume', 'Climate'};
base_path = '~/Dropbox/0osuflow/Paper/LDAV14/fig/matlab';


close all

jet_inv = jet;
jet_inv = jet_inv(length(jet_inv):-1:1,:)
% RMSE
for i=1:length(xvalues)
    figure;
    colormap(jet_inv)
    
    yvalues{i} = abs(yvalues{i});
    bar(100-yvalues{i}'*100);
    if i==1
        h = legend(legends{1}, 'Location', 'NorthWest');
        set(h, 'FontSize', 12.0);
          set(h,'box','off') 
    end
    set(gca, 'FontSize', 36.0);
    h = xlabel({' ';' '; xlabels{1}})
    set(h, 'FontSize', 30.0);
    ylabel(ylabels{1})
    title(titles{i})
    xlim([0.5,size(yvalues{i}, 2)+.5])
    ylim([0, 100])
    set(gca,'YTick', 0:10:100);
    set(gca,'YTickLabel', {0, '', 20, '', 40, '', 60, '', 80, '', 100});

%     set(gca, 'XTickMode', 'manual')
%     set(gca,'XTick', 1:length(xvalues{i}))
%     set(gca,'XTickLabel',xvalues{i})
    my_xticklabels(gca, 1:length(xvalues{i}), xvalues{i})
    saveas(gca, sprintf('%s/%s_%s.eps', base_path, 'Test', titles{i}), 'psc2');
    
end


