
% isabel
TEST = 2;
SAMPLING=13;
switch TEST
    case 0
        case_folder = 'diag3dtime';
        base_path = '/data/flow/diag3dtime';
        true_list_file = sprintf('%s/all.list', base_path);
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=8;        
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;

    case 1
        label = 'isabel'
        base_path = '/data/flow/isabel_all';
        true_list_file = sprintf('%s/all.list', base_path);

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=20; % 20 3125 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
    case 2
        label = 'plume'
        base_path = '/data/flow2/plume126_all';
        true_list_file = sprintf('%s/all.list', base_path);

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=25; % 7200 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
end

sample_w = length(seeds_x)
sample_h = length(seeds_y)
sample_d = length(seeds_z)
sample_t = SAMPLING+1
disp(sprintf('press a key... total points=%d', sample_w*sample_h*sample_d*3))
pause
% get fields
if 1
    [true_list, data_w, data_h, data_d, data_t, scaling] = load_list(true_list_file);
    
    vecAll      = zeros(3, sample_w, sample_h, sample_d, sample_t);
    for t=1:sample_t
       vec      = load_vec(true_list{t});

       % trim
       vecAll(:,:,:,:,t)        = vec(:,seeds_x,seeds_y,seeds_z);
    end

end

% compute linterp field
if 1
    disp('Generating interpolations')
    vecAll_linear = zeros(3, sample_w, sample_h, sample_d, sample_t);
    vecAll_linear2 = zeros(3, sample_w, sample_h, sample_d, sample_t);
    vecAll_bezier = zeros(3, sample_w, sample_h, sample_d, sample_t);
    t_linear  = [1 sample_t];
    t_linear2  = [1 ceil(sample_t/2) sample_t];
    t_ary = 1:sample_t;
    % linterp
    for z=1:sample_d
        for x=1:sample_w
            for y=1:sample_h
                v = squeeze(vecAll(:,x,y,z,:));
                vecAll_linear(:,x,y,z,:) = [interp1q(t_linear', v(1,t_linear)', t_ary')'
                                            interp1q(t_linear', v(2,t_linear)', t_ary')'
                                            interp1q(t_linear', v(3,t_linear)', t_ary')'];     
                vecAll_linear2(:,x,y,z,:) = [interp1q(t_linear2', v(1,t_linear2)', t_ary')'
                                            interp1q(t_linear2', v(2,t_linear2)', t_ary')'
                                            interp1q(t_linear2', v(3,t_linear2)', t_ary')'];     
                vecAll_bezier(:,x,y,z,:) = bezierfit1((t_ary-1), v, (t_ary-1))   ;
            end
        end
    end
% %     err_linear = vecAll_linear - vecAll;
% %     err_linear2 = vecAll_linear2 - vecAll;
    if TEST==1
        vecAll = vecAll(1:2, :,:,:,:);
        vecAll_linear = vecAll_linear(1:2, :,:,:,:);
        vecAll_linear2 = vecAll_linear2(1:2, :,:,:,:);
        vecAll_bezier = vecAll_bezier(1:2, :,:,:,:);
    end
    err_mean_linear = mean(mean(mean(mean(mean(abs(vecAll_linear - vecAll))))))
    err_mean_linear2 = mean(mean(mean(mean(mean(abs(vecAll_linear2 - vecAll))))))
    err_mean_bezier = mean(mean(mean(mean(mean(abs(vecAll_bezier - vecAll))))))

    disp('Running fit test')
    [ks_linear, chi_linear] = normality_test1(vecAll_linear - vecAll)
    [ks_linear2, chi_linear2] = normality_test1(vecAll_linear2 - vecAll)
    [ks_bezier, chi_bezier] = normality_test1(vecAll_bezier - vecAll)
end