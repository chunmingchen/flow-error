function out = gen_seed(filename, xgrid, ygrid, zgrid)

[x,y,z] = meshgrid(xgrid, ygrid, zgrid);
n=size(x,1)*size(x,2)*size(x,3)
x = reshape(x,n,1);
y = reshape(y,n,1);
z = reshape(z,n,1);
% [x y z]
fp=fopen(filename, 'wt');
fprintf(fp, '%d\n', n);
fprintf(fp, '%d %d %d 0\n', [x';y';z']);
fclose(fp);

out = [x y z];

end