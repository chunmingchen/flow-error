function dump_scalar(field, fname)

size(field)
fp=fopen(fname, 'wb');
fwrite(fp, permute(field, [2 1 3]), 'float32');
%fwrite(fp, field, 'float32');
fclose(fp);

end