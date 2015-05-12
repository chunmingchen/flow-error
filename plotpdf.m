base_path = '~/Dropbox/0osuflow/Paper/LDAV14/fig/matlab';

m1=5; m2=10;
s1=1; s2=2;
d1 = makedist('Normal', 'mu', m1, 'sigma', s1);
d2 = makedist('Normal', 'mu', m2, 'sigma', s2);
s3=1/(1/s1+1/s2)
m3=(m1/s1+m2/s2)*s3
d3 = makedist('Normal', 'mu', m3, 'sigma', s3);
x=0:.1:16;
y1=pdf(d1,x);
y2=pdf(d2,x);
y3=pdf(d3,x);

plot( 0:1,0,'r', 'LineWidth', 3);
legend('Intersected Distribution', 'Location', 'Northwest')
legend('boxoff')
hold all        
plot(x,y3,'r','LineWidth', 3)
hold on
plot(x,y1,'b--','LineWidth', 3)
plot(x,y2,'b--','LineWidth', 3)
xlim([1,14])
ylim([0,0.7])
xlabel('X')
ylabel('Probability')
set(gca, 'FontSize', 30.0);
set(gca, 'XTick', []);
hold off
saveas(gca, sprintf('%s/dist.eps', base_path), 'psc2');
