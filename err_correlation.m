err = vecAll_bezier - vecAll;

rmserr =  squeeze(sqrt(mean(err.^2, 5)));

skewerr = squeeze(mean(err.^3, 5));

pfield1 = pfield(1,:,:,:);
rmserr1 = rmserr(1,:,:,:);

% fp = fopen('isabelp_12/kstest_field.raw', 'rb')
% pvalue = fread(fp, 'float32');
% fclose(fp);
 dim = size(pfield1)
% dump_scalar(err, sprintf('kstest_field_%dx%dx%d.raw', dim(1), dim(2), dim(3)));

n = dim(1)*dim(2)*dim(3)*dim(4);
plot(reshape(rmserr1, n, 1), reshape(pfield1, n, 1), '.')
xlabel('stdard error')
ylabel('pvalue')




figure
plot(reshape(skewerr, 200000, 1), reshape(pfield, 200000, 1), '.')
xlabel('skew error')
ylabel('pvalue')
