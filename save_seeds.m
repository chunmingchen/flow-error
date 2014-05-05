% seeds: (4,n)
function save_seeds(filename, seeds)

n = length(seeds);
% [x y z]
fp=fopen(filename, 'wt');
fprintf(fp, '%d\n', n);
fprintf(fp, '%f %f %f %f\n', seeds);
fclose(fp);
end