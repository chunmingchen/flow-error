setenv('PATH', '/home/chenchu/Project/flowvis/tools:/home/chenchu/Project/flowvis/4D/tools:/home/chenchu/bin:/usr/lib/lightdm/lightdm:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games')

skipping = 12;
test = 31;
RUN_TRACING = 1;
switch test
    case 0
        case_folder = 'diag3dtime';
        base_path = '/data/flow/diag3dtime';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=64;
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale=1;
    case 1
        label = 'isabel'
        base_path = '/data/flow/isabel_all';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=20; % 20 3125 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale=1;
    case 11
        label = 'isabel_test'
        base_path = '/data/flow/isabel_all';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

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
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=100; 
        seeds_x = 200:20:250;
        seeds_y = 100:20:150;
        seeds_z = 50;
        scale=1;
    case 2
        label = 'plume'
        base_path = '/data/flow2/plume126_all';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=13; % 7200 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale=1;
    case 21
        
        label = 'plume_vis'
        base_path = '/data/flow2/plume126_all';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=2.5; % 7200 seeds
        seeds_x = 53:SEEDING_STEP:73;
        seeds_y = seeds_x;
        seeds_z = 250;
        scale = 1;
    case 3
        
        label = 'climate'
        base_path = '/data/flow2/smhagos/curvilinear';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=25 % 7200 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        SEEDING_STEP=12; % 5184 * 3 = 15552 seeds
        seeds_z = SEEDING_STEP/2:SEEDING_STEP:D;
        scale = 1;
    case 31
        
        label = 'climate_vis'
        base_path = '/data/flow2/smhagos/curvilinear';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));

        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=50 % 7200 seeds
        seeds_x = SEEDING_STEP/2:SEEDING_STEP:W;
        seeds_y = SEEDING_STEP/2:SEEDING_STEP:H;
        seeds_z = D/2
        scale = 1;
    case 4
        label = 'karman'
        base_path = '/data/flow2/karman/vec/3D/';
        true_list_file = sprintf('%s/all.list', base_path);
        fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);
        dbr_list_file = sprintf('%s/all_%d.list', base_path, ceil(skipping/2));
        
        [files, W, H, D, T, scaling] = load_list(true_list_file);
        SEEDING_STEP=5; % 5184 * 3 = 15552 seeds
        seeds_x = floor(SEEDING_STEP/2):SEEDING_STEP:W-1;
        seeds_y = H/4+floor(SEEDING_STEP/2):SEEDING_STEP:H-H/4;
        seeds_z = 2.5;
        scale = 5;
end
STEP_SIZE = 1/skipping/8
case_folder = sprintf('%s_%d', label, skipping);

% auto
mkdir (case_folder)

% gen seeds
seed_file = sprintf('%s/seeds.txt', case_folder);
gen_seed(seed_file, seeds_x, seeds_y, seeds_z, 0);
% gen_seed(seed_file, 250:20:300, 150:20:200:SEEDING_STEP:H, 50);

%% run 
% truth
traces_true_file = sprintf('%s/true.out', case_folder)
if RUN_TRACING & 1
    cmd = sprintf('parallelPathline_omp %s -loader=1 -saver=2 -limit=3 -omp_nproc=7 -stepsize=.125 -seedfile=%s -out=%s -scale=%f' , ...
        true_list_file, seed_file, traces_true_file, scale)
    system(cmd);
    system(sprintf('convertPathlineVTK %s', traces_true_file));
end

%% sample output
% extract output
traces_true = load_pathlines(traces_true_file);

traces_true_seg = cell(1,10000);
nsegs=0;
for trace = traces_true
    trace1 = trace{1};
    trace1(4,:) = trace1(4,:)/skipping; % N timesteps -> 1 step
    true_t_ary = trace1(4,:);
    for i=true_t_ary(1):true_t_ary(end)-1
        t_ary = i: STEP_SIZE : i+1;
        out_trace = [interp1q(true_t_ary', trace1(1,:)', t_ary')'
                                  interp1q(true_t_ary', trace1(2,:)', t_ary')'
                                  interp1q(true_t_ary', trace1(3,:)', t_ary')'
                                  t_ary];     
                              
        if all( out_trace(1:3,end) >= [0;0;0]) && all( out_trace(1:3,end) <  [W;H;D]) 
            nsegs=nsegs+1;
            traces_true_seg{nsegs} = out_trace;
        end
    end
end
traces_true_seg = traces_true_seg(1:nsegs);

% extrace start and end points
segments = zeros(4, 2, nsegs);
count = 1;
for trace = traces_true_seg
    trace1 = trace{1};
    segments(:, :, count) = [trace1(:,1) trace1(:,end)];
    count=count+1;
end

seeds_fw_file = sprintf('%s/seeds.seg.fw.txt', case_folder)
seeds_bw_file = sprintf('%s/seeds.seg.bw.txt', case_folder)
seeds_dbr_fw_file = sprintf('%s/seeds.seg.dbr.fw.txt', case_folder)
save_seeds(seeds_fw_file, squeeze(segments(:,1,:)));
save_seeds(seeds_bw_file, squeeze(segments(:,2,:)));
segments_dbr = segments;  segments_dbr(4,:,:) = segments_dbr(4,:,:)*2;
save_seeds(seeds_dbr_fw_file, squeeze(segments_dbr(:,1,:)));

% regular forward
traces_fw_rk4_file = sprintf('%s/pathlines.fw.rk4.out', case_folder); 
if RUN_TRACING 
    delete(traces_fw_rk4_file);
    cmd = sprintf('parallelPathline_omp %s -loader=bezier -saver=2 -limit=3 -omp_nproc=7 -stepsize=%.7f -seedfile=%s -out=%s -steps=1 -scale=%f' , ...
        fitted_list_file, STEP_SIZE, seeds_fw_file, traces_fw_rk4_file, skipping*scale)
    system(cmd);
end

% regular forward double sampling rate
traces_dbr_fw_rk4_file = sprintf('%s/pathlines.dbr.fw.rk4.out', case_folder); 
if RUN_TRACING & 1
    delete(traces_dbr_fw_rk4_file);
    cmd = sprintf('parallelPathline_omp %s -loader=1 -saver=2 -limit=3 -omp_nproc=7 -stepsize=%.7f -seedfile=%s -out=%s -steps=2 -scale=%.7f' , ...
        dbr_list_file, 2/skipping/8, seeds_dbr_fw_file, traces_dbr_fw_rk4_file, skipping/2*scale)
    system(cmd);
end

% for demo purpose
out_file = sprintf('%s/pathlines.rk4.out', case_folder); 
if RUN_TRACING & 0
    delete(out_file);
    cmd = sprintf('parallelPathline_omp %s -loader=bezier -saver=2 -limit=3 -omp_nproc=7 -stepsize=%.7f -seedfile=%s -out=%s -scale=%f' , ...
        fitted_list_file, STEP_SIZE, seed_file, out_file, skipping*scale)
    system(cmd);
    system(sprintf('convertPathlineVTK %s', out_file));
end

% regular stochastic
out_file = sprintf('%s/pathlines.ufw.out', case_folder)
if RUN_TRACING & 0
    delete(out_file)
    cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -saver=2 -limit=3 -omp_nproc=7 -stepsize=%.7f -seedfile=%s -out=%s -seedclone=100 -steps=1 -scale=%f' , ...
        fitted_list_file, STEP_SIZE, seeds_fw_file, out_file, skipping*scale)
    system(cmd);
    system(sprintf('convertPathlineVTK %s', out_file));
end


% seeds fw
traces_fw_file = sprintf('%s/pathlines.seg.fw.out', case_folder);  
if RUN_TRACING 
    delete (traces_fw_file);
    cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -gauss_model -saver=2 -limit=3 -omp_nproc=7 -stepsize=%.7f -seedfile=%s -out=%s -steps=1 -scale=%f' , ...
        fitted_list_file, STEP_SIZE, seeds_fw_file, traces_fw_file, skipping*scale)
    system(cmd);
    system(sprintf('convertPathlineVTK %s', traces_fw_file));
end

% seeds bw
traces_bw_file = sprintf('%s/pathlines.seg.bw.out', case_folder); 
if RUN_TRACING 
    delete (traces_bw_file);
    cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -gauss_model -saver=2 -limit=3 -omp_nproc=7 -stepsize=-%.7f -seedfile=%s -out=%s -steps=1 -scale=%f' , ...
        fitted_list_file, STEP_SIZE, seeds_bw_file, traces_bw_file, skipping*scale)
    system(cmd);
    system(sprintf('convertPathlineVTK %s', traces_bw_file));
end

% merge
traces_merged_file = sprintf('%s/pathlines.merged.out', case_folder); 
if RUN_TRACING 
    delete (traces_merged_file);
    cmd = strcat(['build/merge_error  ' , traces_fw_file, ' ', traces_bw_file, ' ', traces_merged_file]);
    system(cmd);
    system(sprintf('convertPathlineVTK %s', traces_merged_file));
end

% statistics
% compare: traces_true_seg vs: linterp, traces_regular,  traces_merged
[traces_merged, traces_merged_std] = load_pathlines(traces_merged_file);
traces_fw_rk4 = load_pathlines(traces_fw_rk4_file);
traces_dbr_fw_rk4 = load_pathlines(traces_dbr_fw_rk4_file);
ntraces = length(traces_true_seg);
% out
off_linterp_list = zeros(skipping*8+1, ntraces);
off_linterp_list2 = zeros(skipping*8+1, ntraces);
off_fw_rk4_list = zeros(skipping*8+1, ntraces);
off_dbr_fw_rk4_list = zeros(skipping*8+1, ntraces);
off_bezier_list = zeros(skipping*8+1, ntraces);
off_merged_list = zeros(skipping*8+1, ntraces);
off_merged_err_list = zeros(3, skipping*8+1, ntraces);
off_merged_std_list = zeros(3, skipping*8+1, ntraces);
count=0;
for i=1:ntraces
    trace_true = traces_true_seg{i};
    t_ary = trace_true(4,:);
    n= length(t_ary);
    
    % linterp
    t2 = [trace_true(4,1) trace_true(4,end)];
    trace_linterp = [interp1q(t2', trace_true(1,[1 n])', t_ary')'
                     interp1q(t2', trace_true(2,[1 n])', t_ary')'
                     interp1q(t2', trace_true(3,[1 n])', t_ary')'
                     t_ary];     
    dist = sqrt(sum((trace_linterp - trace_true).^2));
    off_linterp_list(:,count+1) = dist';
    
    % linterp2
    t2 = trace_true(4,[1 ceil(n/2) n]);
    trace_linterp2 = [interp1q(t2', trace_true(1,[1 ceil(n/2) n])', t_ary')'
                     interp1q(t2', trace_true(2,[1 ceil(n/2) n])', t_ary')'
                     interp1q(t2', trace_true(3,[1 ceil(n/2) n])', t_ary')'
                     t_ary];     
    dist = sqrt(sum((trace_linterp2 - trace_true).^2));
    off_linterp_list2(:,count+1) = dist';
    
    % regular
    trace_fw_rk4 = traces_fw_rk4{i};
    if length(trace_fw_rk4) < n
        continue;
    end
    trace_fw_rk4 = trace_fw_rk4(:,1:n);
    assert(length(trace_fw_rk4) == n);
    
    dist = sqrt(sum((trace_fw_rk4 - trace_true).^2));
    off_fw_rk4_list(:,count+1) = dist';
    
    % regular dbr
    trace_dbr_fw_rk4 = traces_dbr_fw_rk4{i};
    if length(trace_dbr_fw_rk4) < n
        continue;
    end
    trace_dbr_fw_rk4 = trace_dbr_fw_rk4(:,1:n);
    assert(length(trace_dbr_fw_rk4) == n);
    
    dist = sqrt(sum((trace_dbr_fw_rk4 - trace_true).^2));
    off_dbr_fw_rk4_list(:,count+1) = dist';
    
    % merged
    trace_merged = traces_merged{i};
    if isempty(trace_merged) 
        continue;
    end        
    assert(length(trace_merged) == n);
    
    dist = sqrt(sum((trace_merged - trace_true).^2));
    off_merged_list(:,count+1) = dist';
    
    % added 5/22
    % BSpline pathlines
    t2 = (1:length(t_ary))-1;
    trace_bezier = [bezierfit1(t2, trace_true(1,:), t2)
                    bezierfit1(t2, trace_true(2,:), t2)
                    bezierfit1(t2, trace_true(3,:), t2)
                    t_ary];
    dist = sqrt(sum((trace_bezier - trace_true).^2));
    off_bezier_list(:,count+1) = dist';
    
    % merged err/std
    trace_merged_std = traces_merged_std{i};
    assert(length(trace_merged_std) == n);
    
    off_merged_err_list(:,:,count+1) = (trace_merged(1:3,:) - trace_true(1:3,:));
    off_merged_std_list(:,:,count+1) = trace_merged_std;
    
    
    count=count+1;
    
%     plot_traces({trace_true, trace_dbr_fw_rk4, trace_fw_rk4, trace_bezier, trace_merged})
%     disp('pause...')
%     pause
end
off_linterp_list = off_linterp_list(:,1:count);
off_linterp_list2 = off_linterp_list2(:,1:count);
off_fw_rk4_list = off_fw_rk4_list(:,1:count);
off_dbr_fw_rk4_list = off_dbr_fw_rk4_list(:,1:count);
off_bezier_list = off_bezier_list(:,1:count);
off_merged_list = off_merged_list(:,1:count);
off_merged_err_list = off_merged_err_list(:,1:count);
off_merged_std_list = off_merged_std_list(:,1:count);


mean_off_linterp = mean(mean(off_linterp_list)) 
mean_off_linterp2 = mean(mean(off_linterp_list2)) 
mean_off_fw_rk4 = mean(mean(off_fw_rk4_list))
mean_off_dbr_fw_rk4 = mean(mean(off_dbr_fw_rk4_list))
mean_off_merged = mean(mean(off_merged_list ))
mean_off_bezier_list = mean(mean(off_bezier_list))
disp('[mean_off_linterp; mean_off_linterp2;  mean_off_fw_rk4; mean_off_dbr_fw_rk4, mean_off_bezier_list; mean_off_merged]')
[mean_off_linterp; mean_off_linterp2;  mean_off_fw_rk4; mean_off_dbr_fw_rk4; mean_off_bezier_list; mean_off_merged]
n=3*count
s=8
plot(off_merged_std_list(1,:), abs(off_merged_err_list(1,:)),  '.')
hist(off_merged_err_list(1,:) ./ off_merged_std_list(1,:) )
one_times_std = [sum(off_merged_err_list(1,:) <= (off_merged_std_list(1,:)*s)  ) / count
                 sum(off_merged_err_list(2,:) <= (off_merged_std_list(2,:)*s)  ) / count
                 sum(off_merged_err_list(3,:) <= (off_merged_std_list(3,:)*s)  ) / count]
two_times_std = [sum(off_merged_err_list(1,:) <= (off_merged_std_list(1,:)*s*2)  ) / count
                 sum(off_merged_err_list(2,:) <= (off_merged_std_list(2,:)*s*2)  ) / count
                 sum(off_merged_err_list(3,:) <= (off_merged_std_list(3,:)*s*2)  ) / count]
three_times_std = [sum(off_merged_err_list(1,:) <= (off_merged_std_list(1,:)*s*3)  ) / count
                   sum(off_merged_err_list(2,:) <= (off_merged_std_list(2,:)*s*3)  ) / count
                   sum(off_merged_err_list(3,:) <= (off_merged_std_list(3,:)*s*3)  ) / count]
confidence95 = [sum(off_merged_err_list(1,:) <= (off_merged_std_list(1,:)*s*1.96)  ) / count
                   sum(off_merged_err_list(2,:) <= (off_merged_std_list(2,:)*s*1.96)  ) / count
                   sum(off_merged_err_list(3,:) <= (off_merged_std_list(3,:)*s*1.96)  ) / count]      
min_confidence95 = min(confidence95)               
[ntraces count]

% vec_mag(off_merged_err_list) <= (vec_mag(off_merged_std_list)*s) ) / count
%  sum( vec_mag(off_merged_err_list) <= (vec_mag(off_merged_std_list)*2*s) ) / count
%  sum( vec_mag(off_merged_err_list) <= (vec_mag(off_merged_std_list)*3*s) ) / count ]


