
true_list = '/data/flow/isabel_all/2d/all8.list';
fitted_list = '/data/flow/isabel_all/2d/fitted_quad_endpoint/all8.list';
stderr_file = '/data/flow/isabel_all/2d/fitted_quad_endpoint/sampling8_00_stderr.raw';
trace_true = '/data/flow/isabel_all/2d/field_lines_full.txt';
trace_fitted = '/data/flow/isabel_all/2d/fitted_quad_endpoint/field_lines_sampling8.txt';
global data_w data_d data_t data_h 

fitted_list = load_list(fitted_list);

data_t = 9; % sampling+1
t_ary = 0:0.25:data_t-1.25; t_ary = t_ary';

ftrace_true = fopen(trace_true, 'rt');
ftrace_fitted = fopen(trace_fitted, 'rt');
meanerr_fitted_list =[];
meanerr_filtered_list =[];
meanerr_traced_list =[];

while ~feof(ftrace_true)
    % read line
    line = fgets(ftrace_true);
    true_raw_ary = sscanf(line, '%f %f %f %f,', [4 inf]);
    true_raw_ary = true_raw_ary';    
    if true_raw_ary(1,1)<50 || true_raw_ary(1,2)<50 || true_raw_ary(1,1)>450  || true_raw_ary(1,2)>450
        line = fgets(ftrace_fitted);
        continue;
    end

    % even time
    uneven_t_ary = true_raw_ary(:,4);
    true_x_ary = [interp1q(uneven_t_ary, true_raw_ary(:,1), t_ary) ...
                  interp1q(uneven_t_ary, true_raw_ary(:,2), t_ary)];
    true_x_ary = true_x_ary';

    % fit
    [errstd1, errsum1, yy1, a, b, c] = quadfit_endpoint(t_ary,true_x_ary(1,:)',t_ary);
    [errstd2, errsum2, yy2, a, b, c] = quadfit_endpoint(t_ary,true_x_ary(2,:)',t_ary);

    z_ary = [yy1 yy2];
    z_ary=z_ary';
    z_stderr = [errstd1+errsum1; errstd2+errsum2];

    % read line
    line = fgets(ftrace_fitted);
    traced_fitted_x_raw_ary = sscanf(line, '%f %f %f %f,', [4 inf]);
    traced_fitted_x_raw_ary = traced_fitted_x_raw_ary';
    % even time
    traced_fitted_x_ary = [interp1q(traced_fitted_x_raw_ary(:,4), traced_fitted_x_raw_ary(:,1), t_ary) ...
                           interp1q(traced_fitted_x_raw_ary(:,4), traced_fitted_x_raw_ary(:,2), t_ary)];
    traced_fitted_x_ary=traced_fitted_x_ary'; 
    if true_x_ary(:,1)~=traced_fitted_x_ary(:,1)
        error('starting points mismatch')
    end

    % run
    trace_filtering
%     pause
    meanerr_fitted_list(end+1)=meanerr_fitted;
    meanerr_filtered_list(end+1)=meanerr_filtered;
    meanerr_traced_list(end+1)=meanerr_traced;
    
    disp(sprintf('fitted err: %f \nfiltered err: %f \ntraced err: %f\n', ...
        mean(meanerr_fitted_list), mean(meanerr_filtered_list), mean(meanerr_traced_list)))
end

