true_list_file = '/data/flow/isabel_all/2d/all8.list';
fitted_list_file = '/data/flow/isabel_all/2d/fitted_quad_endpoint8/all8.list';
stderr_file = '/data/flow/isabel_all/2d/fitted_quad_endpoint8/sampling8_00_stderr.raw';
global data_w data_d data_t data_h 


% pathline_loader
% output: trace_list

% get fields
if 0
    fitted_list = load_list(fitted_list_file);
    
    vecAll_fit = zeros(3, data_w, data_h, 1, data_t);
    z=2
    for t=1:data_t
       vec_fit = load_vec(fitted_list{t});

       % trim
       vecAll_fit(:,:,:,:,t) = vec_fit(:,:,:,z);
    end

    % get error field => err_fitted
    fin = fopen(stderr_file, 'rb');
    err_fitted = fread(fin, [3 data_w*data_h], 'float32');
    err_fitted = reshape(err_fitted, [3 data_w data_h]);
    fclose(fin);
end

data_t = 9;
dt = 0.1;

for trace = trace_list
    true_raw_ary = trace{1};
    if true_raw_ary(4,end) < data_t-1.1
        continue;
    end
   
    % sampling timing
    t_ary = 0:dt:data_t-dt-1; 
    t_ary = t_ary';
    
    
    % process trace_list - even trace:    
    uneven_t_ary = true_raw_ary(4,:)';
    true_x_ary = [interp1q(uneven_t_ary, true_raw_ary(1,:)', t_ary) ...
                  interp1q(uneven_t_ary, true_raw_ary(2,:)', t_ary)];
    true_x_ary = true_x_ary';

    % fit
    [errstd1, errsum1, yy1, a, b, c, rmserr1] = quadfit_endpoint(t_ary,true_x_ary(1,:)',t_ary);
    [errstd2, errsum2, yy2, a, b, c, rmserr2] = quadfit_endpoint(t_ary,true_x_ary(2,:)',t_ary);

    z_ary = [yy1 yy2];
    z_ary=z_ary';
    z_stderr = zeros(size(z_ary));
    z_stderr(1,:) = rmserr1; %errstd1+errsum1; 
    z_stderr(2,:) = rmserr2; %errstd2+errsum2;

    %%% trace particle
    traced_fitted_x_ary=[true_raw_ary(1,1); true_raw_ary(2,1)];
    pre_x = traced_fitted_x_ary(:,1);
    pre_t = 0;
    pre_x
    for t = t_ary'
        if t==0 
            continue;
        end
        x = runge_kutta4(vecAll_fit, pre_x, pre_t, t-pre_t);
        traced_fitted_x_ary(:, end+1) = x;
        pe_x = x;
        pe_t = t;
        x
    end
    
%     traced_fitted_x_raw_ary=true_raw_ary'; %%[11.000000 11.000000 2.000000 0.000000; 15.784504 11.436985 2.000000 1.000000; 17.084600 11.512656 2.000000 1.250000; 18.414061 11.578578 2.000000 1.500000; 19.773125 11.641154 2.000000 1.750000; 21.158314 11.705994 2.000000 2.000000; 22.555178 11.764643 2.000000 2.250000; 23.968336 11.818457 2.000000 2.500000; 25.417698 11.870479 2.000000 2.750000; 26.905905 11.917877 2.000000 3.000000; 28.409637 11.972345 2.000000 3.250000; 29.909431 12.012804 2.000000 3.500000; 31.437849 12.008892 2.000000 3.750000; 33.005356 11.998375 2.000000 4.000000; 34.573181 11.995104 2.000000 4.250000; 36.156193 11.963609 2.000000 4.500000; 37.761539 11.896452 2.000000 4.750000; 39.372112 11.805254 2.000000 5.000000; 40.969849 11.706592 2.000000 5.250000; 42.560917 11.593554 2.000000 5.500000; 44.161194 11.454815 2.000000 5.750000; 45.768669 11.272562 2.000000 6.000000; 47.391045 11.029262 2.000000 6.250000; 49.026974 10.745399 2.000000 6.500000; 50.657955 10.442777 2.000000 6.750000; 52.266922 10.100547 2.000000 7.000000; 53.861134 9.671535 2.000000 7.250000; 55.444309 9.137972 2.000000 7.500000; 57.030910 8.526555 2.000000 7.750000; 58.645233 7.841821 2.000000 8.000000];
%     traced_fitted_x_ary = [interp1q(traced_fitted_x_raw_ary(:,4), traced_fitted_x_raw_ary(:,1), t_ary) ...
%                            interp1q(traced_fitted_x_raw_ary(:,4), traced_fitted_x_raw_ary(:,2), t_ary)];
%     traced_fitted_x_ary=traced_fitted_x_ary'; 

%     while 1
        trace_filtering

        z_ary = post_x_ary;
        z_stderr = post_x_stderr;
%         pause
%     end

 
    
%     plot(true_raw_ary(1,:), true_raw_ary(2,:))
%     pause
end
