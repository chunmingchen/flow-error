% isabel
skipping = 8;
if 0
    case_folder = 'diag3dtime';
    base_path = '/data/flow/diag3dtime';
    true_list_file = sprintf('%s/all.list', base_path);
    fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);

    W=256; H=256; D=256; T=41;
    SEEDING_STEP=64;
else
    case_folder = 'isabel';
    base_path = '/data/flow/isabel_all';
    true_list_file = sprintf('%s/all.list', base_path);
    fitted_list_file = sprintf('%s/fitted/fitted%d/all_bezier_rms.list', base_path, skipping);

    W=500; H=500; D=100; T=48;
    SEEDING_STEP=100;
end

% auto
mkdir (case_folder)

% gen seeds
seed_file = sprintf('%s/seeds.txt', case_folder);
gen_seed(seed_file, SEEDING_STEP/2:SEEDING_STEP:W, SEEDING_STEP/2:SEEDING_STEP:H, SEEDING_STEP/2:SEEDING_STEP:D);

%% run 
% truth
trace_true_file = sprintf('%s/true.out', case_folder)
if 0
    cmd = sprintf('parallelPathline_omp %s -loader=1 -saver=2 -limit=3 -omp_nproc=7 -stepsize=.125 -seedfile=%s -out=%s' , ...
        true_list_file, seed_file, trace_true_file)
    system(cmd);
    system(sprintf('convertPathlineVTK %s', trace_true_file));
end

%% sample output
% extract output
trace_true = load_pathlines(trace_true_file);

traces_skipped = cell(size(trace_true));
i=1;
nsegs = 0;
for trace = trace_true
    trace = trace{1};
    trace(4,:) = trace(4,:)/skipping; % N files -> 1 file
    true_t_ary = trace(4,:);
    t_ary = true_t_ary(1):1: true_t_ary(end)
%     if t_ary(end) ~= true_t_ary(end)
%         t_ary(end+1) = true_t_ary(end);
%     end
    if length(t_ary)<=1
        continue;
    end
    traces_skipped{i} = [interp1q(true_t_ary', trace(1,:)', t_ary')'
                            interp1q(true_t_ary', trace(2,:)', t_ary')'
                            interp1q(true_t_ary', trace(3,:)', t_ary')'
                            t_ary];    
    
%     plot(trace(1,:), trace(2,:), traces_skipped{i}(1,:), traces_skipped{i}(2,:))
    nsegs  = nsegs + length(t_ary)-1;
%     hold on
    i=i+1;
end
% hold off

segments = zeros(4, 2, nsegs);
count = 1;
for trace = traces_skipped
    trace = trace{1};
    if isempty(trace)
        continue;
    end
    for j=1:size(trace, 2)-1
        segments(:, :, count) = [trace(:,j) trace(:,j+1)];
        count=count+1;
    end
end

seeds_fw_file = sprintf('%s/seeds.seg.fw.txt', case_folder)
seeds_bw_file = sprintf('%s/seeds.seg.bw.txt', case_folder)
save_seeds(seeds_fw_file, squeeze(segments(:,1,:)));
save_seeds(seeds_bw_file, squeeze(segments(:,2,:)));

%% 
% regular forward
cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -gauss_model -saver=2 -limit=3 -omp_nproc=7 -stepsize=0.015625 -seedfile=%s -out=%s/skipping%d.out' , ...
    fitted_list_file, seed_file, case_folder, skipping)
system(cmd);

% regular stochastic
out_file = sprintf('%s/skipping%d.out', case_folder, skipping)
cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -saver=2 -limit=3 -omp_nproc=7 -stepsize=0.015625 -seedfile=%s -out=%s/skipping%d.out -seedclone=100' , ...
    fitted_list_file, seed_file, out_file)
system(cmd);
system(sprintf('convertPathlineVTK %s', out_file));


% seeds fw
trace_fw_file = sprintf('%s/pathlines.seg.fw.out', case_folder)
cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -gauss_model -saver=2 -limit=3 -omp_nproc=7 -stepsize=0.015625 -seedfile=%s -out=%s -steps=1' , ...
    fitted_list_file, seeds_fw_file, trace_fw_file)
system(cmd);
system(sprintf('convertPathlineVTK %s', trace_fw_file));

% seeds bw
trace_bw_file = sprintf('%s/pathlines.seg.bw.out', case_folder)
cmd = sprintf('parallelPathline_omp %s -loader=bezier_random -gauss_model -saver=2 -limit=3 -omp_nproc=7 -stepsize=-0.015625 -seedfile=%s -out=%s -steps=1' , ...
    fitted_list_file, seeds_bw_file, trace_bw_file)
system(cmd);
system(sprintf('convertPathlineVTK %s', trace_bw_file));

% merge
trace_merged_file = sprintf('%s/pathlines.merged.out', case_folder);
cmd = strcat(['build/merge_error  ' , trace_fw_file, ' ', trace_bw_file, ' ', trace_merged_file]);
system(cmd);
system(sprintf('convertPathlineVTK %s', trace_merged_file));

% statistics
trace_merged = load_pathlines(trace_merged_file);

