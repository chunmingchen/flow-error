function save_vec(filename, vec)

disp(sprintf('saving %s...', filename))
fout=fopen(filename, 'wb')
if (fout < 0)
    error('File not found');
end
s = size(vec);

dim = [s(2) s(3) s(4)]; % first dim of s is 3
dim

fwrite(fout, dim, 'int32');
fwrite(fout, vec, 'float32');

fclose(fout);

end