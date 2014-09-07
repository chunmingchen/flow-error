setenv('PATH', '/home/chenchu/Project/flowvis/tools:/home/chenchu/Project/flowvis/4D/tools:/home/chenchu/bin:/usr/lib/lightdm/lightdm:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games')

sampling = 1;
test = 1;

RUN_TRACING = 1;
switch test
    case 0
        case_folder = 'diag3dtime';
        base_path = '/data/flow/diag3dtime';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=64;
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale=1;
    case 1
        label = 'isabel'
        base_path = '/data/flow/isabel_all';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=1; % 20 3125 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale=1;
    case 11
        label = 'isabel_test'
        base_path = '/data/flow/isabel_all';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=20; 
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = 50;
        scale=1;
        %seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
    case 12
        label = 'isabel_vis'
        base_path = '/data/flow/isabel_all';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=100; 
        seeds_x = 200:20:250;
        seeds_y = 100:20:150;
        seeds_z = 50;
        scale=1;
    case 2
        label = 'plume'
        base_path = '/data/flow2/plume252_all';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=1; % 7200 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale=1;
    case 21
        
        label = 'plume_vis'
        base_path = '/data/flow2/plume126_all';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=2.5; % 7200 seeds
        seeds_x = 53:SEEDING_STEP:73;
        seeds_y = seeds_x;
        seeds_z = 250;
        scale = 1;
    case 3
        
        label = 'climate'
        base_path = '/data/flow2/smhagos/curvilinear';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=1 % 7200 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale = 1;
    case 31
        
        label = 'climate_vis'
        base_path = '/data/flow2/smhagos/curvilinear';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=50 % 7200 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = D/2
        scale = 1;
    case 4
        label = 'karman'
        base_path = '/data/flow2/karman/vec/3D/';
        true_list_file = sprintf('%s/all2.list', base_path);
        %fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, sampling);
        %dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(sampling/2));
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=5; % 5184 * 3 = 15552 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W-1;
        seeds_y = H/4+floor(SEEDING_STEP/2):SEEDING_STEP:H-H/4;
        seeds_z = 2.5;
        scale = 5;
end
%STEP_SIZE = 1/sampling/8
%case_folder = sprintf('%s_%d', label, sampling);
case_folder = sprintf('%s', label);

% auto
mkdir (case_folder)

% gen seeds
seed_file = sprintf('%s/seeds.txt', case_folder);
if 0
    gen_seed(seed_file, seeds_x, seeds_y, seeds_z, 0, false);
    % gen_seed(seed_file, 250:20:300, 150:20:200:SEEDING_STEP:H, 50);
end

%% run 
% truth
traces_true_file = sprintf('%s/trace_all.out', case_folder)
omp_nproc=4
if RUN_TRACING & 1
    cmd = sprintf('parallelPathline_omp %s -loader=1 -saver=0 -limit=3 -omp_nproc=%d -stepsize=.25 -seedfile=%s -out=%s -scale=%f -maxT=1 -nproc=1' , ...
        true_list_file, omp_nproc, seed_file, traces_true_file, scale)

    system(cmd);
    %system(sprintf('convertPathlineVTK %s', traces_true_file));
end

