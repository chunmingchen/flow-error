test=3

switch (test)
    case 1
        err = vecAll_bezier - vecAll;

        rmserr =  squeeze(sqrt(mean(err.^2, 5)));
        pfield = reshape(pfield, [100,100,20]);
    case 2
        % 
        % skewerr = squeeze(mean(err.^3, 5));

        % fp = fopen('isabelp_12/kstest_field.raw', 'rb')
        % pvalue = fread(fp, 'float32');
        % fclose(fp);
        pfield = reshape(pfield, [100,100,20]);

    case 3
        err = vecAll_fitn1 - vecAll;
        rmserr =  squeeze(sqrt(mean(err.^2, 5)));
        pfield = pfield1
        
end

s=size(rmserr);
n = s(1)*s(2)*s(3)*s(4);

plot(reshape(rmserr, n, 1), reshape(pfield, n, 1), '.')
xlabel('stdard error')
ylabel('pvalue')



% 
% figure
% plot(reshape(skewerr, 200000, 1), reshape(pfield, 200000, 1), '.')
% xlabel('skew error')
% ylabel('pvalue')
