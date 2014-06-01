function out = gen_seed(filename, xgrid, ygrid, zgrid, tgrid)

[x,y,z] = meshgrid(xgrid, ygrid, zgrid);
n=size(x,1)*size(x,2)*size(x,3)
x = reshape(x,n,1);
y = reshape(y,n,1);
z = reshape(z,n,1);
% [x y z]
fp=fopen(filename, 'wt');
fprintf(fp, '%d\n', n);
out=[];
for t=tgrid
    tary = ones(n,1)*t;
    fprintf(fp, '%g %g %g %g\n', [x';y';z';tary']);
    out = [out; [x y z tary]];
end
fclose(fp);


end