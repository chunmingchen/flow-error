function save_list(path, files, w, h, d, t)
fname = sprintf('%s/all.list', path);
fname
fout=fopen(fname, 'wt');
fprintf(fout, '%d %d %d %d \t# w h d t\n', w, h, d, t);
fprintf(fout, '1 \t# scaling\n');
for file=files
    [pathstr,name,ext] = fileparts(file{:});
    fprintf(fout, '%s%s\n', name, ext);
end
fclose(fout);

end