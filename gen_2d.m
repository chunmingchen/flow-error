list = '/data/flow/isabel_all/all.list';
out_list = '/data/flow/isabel_all/2d';
GET_Z = 50
RESIZE_Z = 5

global data_w data_d data_t data_h 
files = load_list(list);

out_files = {};
for t=1:data_t
    file = files{t};
    
    vec = load_vec(file);
    vec2d = vec(:,:,:,GET_Z);
    vec2d(3,:,:,:)=0;
    
    % replicate
    s = size(vec2d); s(4)=RESIZE_Z;
    vec3d = zeros(s);
    for i=1:RESIZE_Z
        vec3d(:,:,:,i) = vec2d;
    end
    
    [pathstr,name,ext] = fileparts(files{t});
    out_file = sprintf('%s/%s_2d%s', out_list, name, ext);
    if (strcmp( out_file, file )==1) 
        error('write to the same file');
        break
    end
    save_vec(out_file, vec3d);
    
    out_files = {out_files{:}, out_file};
end

save_list(out_list, out_files, data_w, data_h, RESIZE_Z, data_t);
