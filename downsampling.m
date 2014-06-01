if 1
        base_path = '/data/flow2/smhagos';
        true_list_file = sprintf('%s/all.list', base_path);
        downsample_step = 25
end    
if 1
    [true_list, data_w, data_h, data_d, data_t, scaling] = load_list(true_list_file);
    
    for t=1:data_t
        system('downsampling %s %d
    end

end