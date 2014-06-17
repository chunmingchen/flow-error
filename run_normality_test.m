
% isabel
TEST = 101
SAMPLING=47;
TUPLES=3; % default
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
        %SEEDING_STEP=16; % 20 3125 * 2 = 9375 seeds %<< paper
        SEEDING_STEP=5; % 20 3125 * 2 = 9375 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
        
    case 2
        label = 'plume'
        base_path = '/data/flow2/plume126_all';
        true_list_file = sprintf('%s/all.list', base_path);

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=13; % 4100 * 3 = 11700 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
    case 3
        label = 'smhagos'
        base_path = '/data/flow2/smhagos/curvilinear';
        true_list_file = sprintf('%s/all.list', base_path);
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=25; % 5184 * 3 = 15552 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        SEEDING_STEP=12; % 5184 * 3 = 15552 seeds
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
    case 4
        label = 'karman'
        base_path = '/data/flow2/karman/vec/3D/';
        true_list_file = sprintf('%s/all.list', base_path);
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=5; % 5184 * 3 = 15552 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = 1;
    case 5 
        label = 'tornado'
        base_path = '/data/flow2/tornado/500_50x';
        true_list_file = sprintf('%s/all.list', base_path);
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=8; % 4096 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
    
        %%%%%%%%%%%%%% scalar
    case 101
        label = 'isabelp'
        base_path = '/data/flow2/isabel_all';
        true_list_file = sprintf('%s/all.list', base_path);

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=5; % 20 3125 * 2 = 9375 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = floor(SEEDING_STEP/2):SEEDING_STEP:H;
        seeds_z = floor(SEEDING_STEP/2):SEEDING_STEP:D;
        TUPLES=1
end

sample_w = length(seeds_x)
sample_h = length(seeds_y)
sample_d = length(seeds_z)
sample_t = SAMPLING+1
disp(sprintf('press a key... total points=%d *%d=%d', sample_w*sample_h*sample_d, TUPLES, sample_w*sample_h*sample_d*TUPLES))
pause
% get fields
if 1
    [true_list, data_w, data_h, data_d, data_t, scaling] = load_list(true_list_file);
    
    vecAll      = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    for t=1:sample_t
        if TUPLES==1
            vec     = load_scalar(true_list{t}, W, H, D);
        else
            vec     = load_vec(true_list{t});
        end

        % trim
        vecAll(:,:,:,:,t)        = vec(:,seeds_x,seeds_y,seeds_z);
    end

end

% compute linterp field
if 1
    disp('Generating interpolations')
    vecAll_linear = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    vecAll_linear2 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    vecAll_bezier = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    t_linear  = [1 sample_t];
    t_linear2  = [1 ceil(sample_t/2) sample_t];
    t_ary = 1:sample_t;
    % linterp
    for z=1:sample_d
        for x=1:sample_w
            for y=1:sample_h
                v = squeeze(vecAll(:,x,y,z,:));
                if TUPLES==3
                    vecAll_linear(:,x,y,z,:) = [interp1q(t_linear', v(1,t_linear)', t_ary')'
                                                interp1q(t_linear', v(2,t_linear)', t_ary')'
                                                interp1q(t_linear', v(3,t_linear)', t_ary')'];     
                    vecAll_linear2(:,x,y,z,:) = [interp1q(t_linear2', v(1,t_linear2)', t_ary')'
                                                interp1q(t_linear2', v(2,t_linear2)', t_ary')'
                                                interp1q(t_linear2', v(3,t_linear2)', t_ary')'];     
                    vecAll_bezier(:,x,y,z,:) = bezierfit1((t_ary-1), v, (t_ary-1))   ;
                else
                    vecAll_linear(:,x,y,z,:) = interp1q(t_linear', v(t_linear), t_ary')';     
                    vecAll_linear2(:,x,y,z,:) = interp1q(t_linear2', v(t_linear2), t_ary')';     
                    vecAll_bezier(:,x,y,z,:) = bezierfit1((t_ary-1), v', (t_ary-1))   ;
                end
            end
        end
    end
% %     err_linear = vecAll_linear - vecAll;
% %     err_linear2 = vecAll_linear2 - vecAll;
    if TEST==1 || TEST==3 || TEST==4 
        vecAll = vecAll(1:2, :,:,:,:);
        vecAll_linear = vecAll_linear(1:2, :,:,:,:);
        vecAll_linear2 = vecAll_linear2(1:2, :,:,:,:);
        vecAll_bezier = vecAll_bezier(1:2, :,:,:,:);
        disp(sprintf('Actual total points=%d', sample_w*sample_h*sample_d*2))
    end
    err_rms_linear = sqrt(mean(mean(mean(mean(mean((vecAll_linear - vecAll).^2))))))
    err_rms_linear2 = sqrt(mean(mean(mean(mean(mean((vecAll_linear2 - vecAll).^2))))))
    err_rms_bezier = sqrt(mean(mean(mean(mean(mean((vecAll_bezier - vecAll).^2))))))
        
    err_mean_linear = mean(mean(mean(mean(mean(mean(vecAll_linear - vecAll))))))
    err_mean_linear2 = mean(mean(mean(mean(mean(mean(vecAll_linear2 - vecAll))))))
    err_mean_bezier = mean(mean(mean(mean(mean(mean(vecAll_bezier - vecAll))))))
    
%     nsr_mean_linear = mean(mean(mean(mean(mean(mean((vecAll_linear - vecAll)./vecAll))))))
%     nsr_mean_linear2 = mean(mean(mean(mean(mean(mean((vecAll_linear2 - vecAll)./vecAll))))))
%     nsr_mean_bezier = mean(mean(mean(mean(mean(mean((vecAll_bezier - vecAll)./vecAll))))))
%     %hist(reshape(mean(vecAll_bezier - vecAll, 5), 11700,1))
%      hist(reshape(mean(vecAll_bezier(1,:,:,:,:) - vecAll(1,:,:,:,:), 5)./mean(vecAll(2,:,:,:,:),5), numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%      hist(reshape(log(abs(mean(vecAll_bezier(1,:,:,:,:) - vecAll(1,:,:,:,:), 5)./mean(vecAll(2,:,:,:,:),5)))/log(10), numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%      hist(reshape(log(abs(mean(vecAll_bezier(2,:,:,:,:) - vecAll(2,:,:,:,:), 5)))/log(10)  , numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%       hist(reshape(log(sqrt(mean((vecAll_bezier(2,:,:,:,:) - vecAll(2,:,:,:,:)).^2, 5)))/log(10)  , numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%     
%     disp('Running fit test')
%     [ks_linear, chi_linear] = normality_test1(vecAll_linear - vecAll)
%     [ks_linear2, chi_linear2] = normality_test1(vecAll_linear2 - vecAll)
%     [ks_bezier, chi_bezier] = normality_test1(vecAll_bezier - vecAll)
    
    [ks_bezier, chi_bezier, conf_ratio, pfield] = normality_test1(vecAll_bezier - vecAll)   ; 
    [ks_bezier, chi_bezier, conf_ratio]
    dump_scalar(squeeze(pfield(1,:,:,:,:)) , sprintf('%s_%d/kstest_field.raw', label, SAMPLING));
    [ks_linear2, chi_linear2] = normality_test1(vecAll_linear2 - vecAll)
    [ks_linear, chi_linear] = normality_test1(vecAll_linear - vecAll)
end
err_rms = [err_rms_linear; err_rms_linear2;  err_rms_bezier]
err_mean = [err_mean_linear; err_mean_linear2; err_mean_bezier]
kstest = [ks_linear;  ks_linear2;    ks_bezier]
chitest = [chi_linear;  chi_linear2;    chi_bezier]
conf_ratio