err = vecAll_bezier - vecAll;

rmserr =  squeeze(sqrt(mean(err.^2, 5)));

skewerr = squeeze(mean(err.^3, 5));

% fp = fopen('isabelp_12/kstest_field.raw', 'rb')
% pvalue = fread(fp, 'float32');
% fclose(fp);

pfield = reshape(pfield, [100,100,20]);


plot(reshape(rmserr, 200000, 1), reshape(pfield, 200000, 1), '.')
xlabel('stdard error')
ylabel('pvalue')




figure
plot(reshape(skewerr, 200000, 1), reshape(pfield, 200000, 1), '.')
xlabel('skew error')
ylabel('pvalue')
