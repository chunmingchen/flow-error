
% isabel
TEST = 3
SAMPLING=100;
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
        SEEDING_STEP=20; % 20 3125 * 2 = 9375 seeds
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
    case 3.1
        label = 'smhagos'
        base_path = '/data/smhagos/franklin/curvelinear';
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
    case 6
        label = 'eddy'
        base_path = '/data/flow/eddy';
        true_list_file = sprintf('%s/all.list', base_path);
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=10; % 4096 seeds
        seeds_x = ceil(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = 1;
        seeds_z = 1;
    case 7
        label = 'turbine'
        base_path = '/data/flow2/turbine_Stg/zDIR.P3D.rel.6201-11001/block0/'
        true_list_file = sprintf('%s/all.list', base_path);
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=100;
        seeds_x = ceil(SEEDING_STEP/2):SEEDING_STEP:W;
        seeds_y = 1;
        seeds_z = 1;
    
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
if 0
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

% These datasets behave more like 2D
if TEST==1 || TEST==3 || TEST==4 || TEST==6
    TUPLES=2 
    vecAll = vecAll(1:2, :,:,:,:);
end

% compute linterp field
if 1
    disp('Generating interpolations')
    vecAll_linear = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    vecAll_linear2 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    vecAll_bezier = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    vecAll_bezier3 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t); % cubic
    vecAll_bezier4 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t); % quartic
    vecAll_bezier5 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t); % quintic
    vecAll_fitn1 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);
    vecAll_fitn2 = zeros(TUPLES, sample_w, sample_h, sample_d, sample_t);    
    t_linear  = [1 sample_t];
    t_linear2  = [1 ceil(sample_t/2) sample_t];
    t_ary = 1:sample_t;
    % linterp
    for z=1:sample_d
        for x=1:sample_w
            for y=1:sample_h
                v = squeeze(vecAll(:,x,y,z,:));
                if TUPLES==1
                    v=v';
                end
                for i=1:TUPLES
                    vecAll_linear(i,x,y,z,:) = interp1q(t_linear', v(i,t_linear)', t_ary')';
                    vecAll_linear2(i,x,y,z,:) = interp1q(t_linear2', v(i,t_linear2)', t_ary')'; 
                    vecAll_bezier(i,x,y,z,:) = bezierfit1((t_ary-1), v(i,:), (t_ary-1))   ;
                    vecAll_bezier3(i,x,y,z,:) = bezierfitcubic(t_ary, v(i,:), t_ary);
                    vecAll_bezier4(i,x,y,z,:) = bezierfitquartic(t_ary, v(i,:), t_ary);
                    vecAll_bezier5(i,x,y,z,:) = bezierfitQuintic(t_ary, v(i,:), t_ary);
                    p = polyfit(t_ary, v(i,:), 1); vecAll_fitn1(i,x,y,z,:) = polyval(p, t_ary);
                    p = polyfit(t_ary, v(i,:), 2); vecAll_fitn2(i,x,y,z,:) = polyval(p, t_ary);                    
                end    
            end
        end
    end
% %     err_linear = vecAll_linear - vecAll;
% %     err_linear2 = vecAll_linear2 - vecAll;
    % this computation may be wrong:
    err_rms_linear = sqrt(mean(mean(mean(mean(mean((vecAll_linear - vecAll).^2))))))
    err_rms_linear2 = sqrt(mean(mean(mean(mean(mean((vecAll_linear2 - vecAll).^2))))))
    err_rms_bezier = sqrt(mean(mean(mean(mean(mean((vecAll_bezier - vecAll).^2))))))
    err_rms_bezier3 = sqrt(mean(mean(mean(mean(mean((vecAll_bezier3 - vecAll).^2))))))
    err_rms_bezier4 = sqrt(mean(mean(mean(mean(mean((vecAll_bezier4 - vecAll).^2))))))
    err_rms_bezier5 = sqrt(mean(mean(mean(mean(mean((vecAll_bezier5 - vecAll).^2))))))
    err_rms_fitn1 = sqrt(mean(mean(mean(mean(mean((vecAll_fitn1 - vecAll).^2))))))
    err_rms_fitn2 = sqrt(mean(mean(mean(mean(mean((vecAll_fitn2 - vecAll).^2))))))
        
    err_mean_linear = mean(mean(mean(mean(mean(mean(vecAll_linear - vecAll))))))
    err_mean_linear2 = mean(mean(mean(mean(mean(mean(vecAll_linear2 - vecAll))))))
    err_mean_bezier = mean(mean(mean(mean(mean(mean(vecAll_bezier - vecAll))))))
    err_mean_bezier3 = mean(mean(mean(mean(mean(mean(vecAll_bezier3 - vecAll))))))
    err_mean_bezier4 = mean(mean(mean(mean(mean(mean(vecAll_bezier4 - vecAll))))))
    err_mean_bezier5 = mean(mean(mean(mean(mean(mean(vecAll_bezier5 - vecAll))))))
    err_mean_fitn1 = mean(mean(mean(mean(mean(mean(vecAll_fitn1 - vecAll))))))
    err_mean_fitn2 = mean(mean(mean(mean(mean(mean(vecAll_fitn2 - vecAll))))))
    
%     nsr_mean_linear = mean(mean(mean(mean(mean(mean((vecAll_linear - vecAll)./vecAll))))))
%     nsr_mean_linear2 = mean(mean(mean(mean(mean(mean((vecAll_linear2 - vecAll)./vecAll))))))
%     nsr_mean_bezier = mean(mean(mean(mean(mean(mean((vecAll_bezier - vecAll)./vecAll))))))
%     %hist(reshape(mean(vecAll_bezier - vecAll, 5), 11700,1))
%      hist(reshape(mean(vecAll_bezier(1,:,:,:,:) - vecAll(1,:,:,:,:), 5)./mean(vecAll(2,:,:,:,:),5), numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%      hist(reshape(log(abs(mean(vecAll_bezier(1,:,:,:,:) - vecAll(1,:,:,:,:), 5)./mean(vecAll(2,:,:,:,:),5)))/log(10), numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%      hist(reshape(log(abs(mean(vecAll_bezier(2,:,:,:,:) - vecAll(2,:,:,:,:), 5)))/log(10)  , numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%       hist(reshape(log(sqrt(mean((vecAll_bezier(2,:,:,:,:) - vecAll(2,:,:,:,:)).^2, 5)))/log(10)  , numel(vecAll)/size(vecAll,5)/size(vecAll,1) ,1),30)
%     
    disp('Running fit test')
    
%     [ks_bezier_search, chi_bezier_search, conf_ratio_bezier_search, pfield_bezier_search, err_rms_bezier_search] = normality_search(vecAll);
    [ks_bezier, chi_bezier, conf_ratio_bezier, pfield_bezier] = normality_test1(vecAll_bezier - vecAll)   ; 
%     [ks_bezier, chi_bezier, conf_ratio_bezier]
    [ks_bezier3, chi_bezier3, conf_ratio_bezier3, pfield_bezier3] = normality_test1(vecAll_bezier3 - vecAll, 2)   ; 
    [ks_bezier4, chi_bezier4, conf_ratio_bezier4, pfield_bezier4] = normality_test1(vecAll_bezier4 - vecAll, 3)   ; 
    %dump_scalar(squeeze(pfield(1,:,:,:,:)) , sprintf('%s_%d/kstest_field.raw', label, SAMPLING));
%     [ks_linear, chi_linear] = normality_test1(vecAll_linear - vecAll)
%     [ks_linear2, chi_linear2] = normality_test1(vecAll_linear2 - vecAll)
%     [ks_fitn1, chi_fitn1, conf_ratio_fitn1, pfield_fitn1] = normality_test1(vecAll_fitn1 - vecAll, 2)
%     [ks_fitn2, chi_fitn2, conf_ratio_fitn2, pfield_fitn2] = normality_test1(vecAll_fitn2 - vecAll, 3)
end
err_rms = [err_rms_linear; err_rms_linear2; err_rms_bezier; err_rms_bezier3; err_rms_bezier4; err_rms_fitn1; err_rms_fitn2]
err_mean = [err_mean_linear; err_mean_linear2; err_mean_bezier; err_mean_bezier3; err_mean_bezier4; err_mean_fitn1; err_mean_fitn2]
kstest = [ks_linear;  ks_linear2;   ks_bezier; ks_bezier3; ks_bezier4;  ks_fitn1;  ks_fitn2]
chitest = [chi_linear;  chi_linear2;  chi_bezier; chi_bezier3; chi_bezier4;  chi_fitn1;  chi_fitn2]
conf_ratio = [conf_ratio_bezier; conf_ratio_bezier3; conf_ratio_bezier4; conf_ratio_fitn1; conf_ratio_fitn2]