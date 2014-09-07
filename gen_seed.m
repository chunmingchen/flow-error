function out = gen_seed(filename, xgrid, ygrid, zgrid, tgrid, bOutput)

[x,y,z] = meshgrid(xgrid, ygrid, zgrid);
n=size(x,1)*size(x,2)*size(x,3)
x = reshape(x,1,n);
y = reshape(y,1,n);
z = reshape(z,1,n);
% [x y z]
fp=fopen(filename, 'wt');
fprintf(fp, '%d\n', n);
out=[];
for t=tgrid
    tary = ones(n,1)*t;
    fprintf(fp, '%g %g %g %g\n', [x;y;z;tary']);
    if bOutput
        out = [out; [x y z tary]];
    end
end
fclose(fp);


end