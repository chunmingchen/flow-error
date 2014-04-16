SKIP=10

load_list
vecs=cell(ceil(data_t/SKIP),1);
for i=1:SKIP:data_t
    for j=i:i+SKIP-1
        vecs{j} = load_vec(files{i});
    end
end
